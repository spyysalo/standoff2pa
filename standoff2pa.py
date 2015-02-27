#!/usr/bin/env python

import sys
import re
import fileinput
import json
import codecs

from collections import defaultdict
from itertools import count
from os import path

DEFAULT_ENCODING='utf-8'
DEFAULT_DB='PubMed'

def argparser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source', metavar='DB', default=DEFAULT_DB,
                        help='Source database (default %s)' % DEFAULT_DB)
    parser.add_argument('-o', '--output', metavar='DIR', default=None,
                        help='Output directory')
    parser.add_argument('-t', '--textdir', metavar='DIR', default=None,
                        help='Text directory')
    parser.add_argument('file', metavar='FILE', nargs='+',
                        help='File(s) to convert')

    return parser

def new_id(prefix, ann_by_id):
    for i in count(1):
        if prefix+str(i) not in ann_by_id:
            return prefix+str(i)

class Annotation(object):
    """Base class for annotations with ID and type."""

    def __init__(self, id_, type_):
        self.id = id_
        self.type = type_

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        raise NotImplementedError

    def verify_text(self, text):
        """Verify reference text for textbound annotations."""
        pass

    def pa_id(self):
        """Return the PubAnnotation ID for the annotation."""
        return self.id

    def to_pubannotation(self, ann_by_id):
        raise NotImplementedError

    STANDOFF_RE = None

    @classmethod
    def from_standoff(cls, line):
        if cls.STANDOFF_RE is None:
            raise NotImplementedError
        m = cls.STANDOFF_RE.match(line)
        if not m:
            raise ValueError('Failed to parse "%s"' % line)
        return cls(*m.groups())

class Textbound(Annotation):
    """Textbound annotation representing entity mention or event trigger."""

    def __init__(self, id_, type_, spans, text):
        super(Textbound, self).__init__(id_, type_)
        self.spans = spans
        self.text = text

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        spans = []
        for span in self.spans.split(';'):
            start, end = span.split(' ')
            spans.append((int(start), int(end)))
        return spans

    def verify_text(self, text):
        offset = 0
        for start, end in self.get_spans():
            endoff = offset + (end-start)
            assert text[start:end] == self.text[offset:endoff], \
                'Error: text mismatch: "%s" vs. "%s"' % \
                (text[start:end], self.text[offset:endoff])
            offset = endoff + 1

    def to_pubannotation(self, ann_by_id):        
        spans = self.get_spans()
        start, end = spans[0][0], spans[-1][1]
        if len(spans) > 1:
            print >> sys.stderr, 'Warning: flattening span %s to %d-%d' % \
                (self.spans, start, end)
        doc = {
            'id': self.pa_id(),
            'obj': self.type,
            'span': { 'begin': start, 'end': end },
        }
        return {
            'denotations': [doc],
        }

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+) (\d+ \d+(?:;\d+ \d+)*)\t(.*)$')

class Relation(Annotation):
    """Typed binary relation annotation."""

    def __init__(self, id_, type_, args):
        super(Relation, self).__init__(id_, type_)
        self.args = args

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        if ann_by_id is None:
            print >> sys.stderr, 'Relation.get_spans: missing ann_by_id'
            return []
        else:
            arg1, arg2 = self.get_args()
            a1, a2 = arg1[1], arg2[1]
            return (ann_by_id[a1].get_spans(ann_by_id) +
                    ann_by_id[a2].get_spans(ann_by_id))

    def get_args(self):
        a1, a2 = self.args.split(' ')
        a1key, a1val = a1.split(':', 1)
        a2key, a2val = a2.split(':', 1)
        return ((a1key, a1val), (a2key, a2val))
    
    def to_pubannotation(self, ann_by_id):
        arg1, arg2 = self.get_args()
        doc = {
            'id': self.pa_id(),
            'pred': self.type,
            'subj': arg1[1],
            'obj': arg2[1],
        }
        return {
            'relations': [doc],
        }

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+) (\S+:\S+ \S+:\S+)$')

class Event(Annotation):
    """Typed, textbound event annotation."""

    def __init__(self, id_, type_, trigger, args):
        super(Event, self).__init__(id_, type_)
        self.trigger = trigger
        self.args = args

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        if ann_by_id is None:
            print >> sys.stderr, 'Event.get_spans: missing ann_by_id'
            return []
        else:
            return ann_by_id[self.trigger].get_spans(ann_by_id)

    def get_args(self):
        return [a.split(':', 1) for a in self.args.split(' ')]        

    def pa_id(self):
        """Return the PubAnnotation ID for the annotation."""
        # Events are represented using their triggers only.
        return self.trigger

    def to_pubannotation(self, ann_by_id):
        relations = []
        for key, val in self.get_args():
            rid = new_id('R', ann_by_id)
            ann_by_id[rid] = None # reserve
            doc = {
                'id': rid,
                'pred': 'has'+key,
                'subj': ann_by_id[val].pa_id(),
                'obj': self.trigger,
            }
            relations.append(doc)
        return {
            'relations': relations,
        }        

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+):(\S+) (\S+:\S+ ?)*$')

class Normalization(Annotation):
    """Reference relating annotation to external resource."""

    def __init__(self, id_, type_, arg, ref, text):
        super(Normalization, self).__init__(id_, type_)
        self.arg = arg
        self.ref = ref
        self.text = text

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        if ann_by_id is None:
            print >> sys.stderr, 'Normalization.get_spans: missing ann_by_id'
            return []
        else:
            return ann_by_id[self.arg].get_spans(ann_by_id)

    def to_pubannotation(self, ann_by_id):
        # map to denotation + relation using Pierre Zweigenbaum's approach:
        # https://github.com/linkedannotation/blah2015/wiki/PubAnnotation-system
        spans = self.get_spans(ann_by_id)
        start, end = spans[0][0], spans[-1][1]
        denotation = {
            'id': self.pa_id(),
            'obj': self.ref,
            'span': { 'begin': start, 'end': end },
        }
        rid = new_id('R', ann_by_id)
        ann_by_id[rid] = None # reserve
        relation = {
            'id': rid,
            'pred': 'Normalization',
            'subj': ann_by_id[self.arg].pa_id(),
            'obj': self.pa_id(),
        }
        return {
            'denotations': [denotation],
            'relations': [relation],
        }

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+) (\S+) (\S+:\S+)\t?(.*)$')

class Attribute(Annotation):
    """Attribute with optional value associated with another annotation."""

    def __init__(self, id_, type_, arg, val):
        super(Attribute, self).__init__(id_, type_)
        self.arg = arg
        self.val = val

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        if ann_by_id is None:
            print >> sys.stderr, 'Attribute.get_spans: missing ann_by_id'
            return []
        else:
            return ann_by_id[self.arg].get_spans(ann_by_id)

    def to_pubannotation(self, ann_by_id):
        pred = self.type + (self.val if self.val is not None else '')
        doc = {
            'id': self.pa_id(),
            'pred': pred,
            'obj': ann_by_id[self.arg].pa_id(),
        }
        return {
            'modifications': [doc],
        }

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+) (\S+) ?(\S*)$')

class Comment(Annotation):
    """Typed free-form text comment associated with another annotation."""

    def __init__(self, id_, type_, arg, text):
        super(Comment, self).__init__(id_, type_)
        self.arg = arg
        self.text = text

    def get_spans(self, ann_by_id=None):
        """Return list of associated (start, end) spans."""
        if ann_by_id is None:
            print >> sys.stderr, 'Comment.get_spans: missing ann_by_id'
        else:
            return ann_by_id[self.arg].get_spans(ann_by_id)

    def to_pubannotation(self, ann_by_id):
        # map to denotation using span associated with target
        spans = self.get_spans(ann_by_id)
        start, end = spans[0][0], spans[-1][1]
        doc = {
            'id': self.pa_id(),
            'obj': self.text,
            'span': { 'begin': start, 'end': end },
        }
        return {
            'denotations': [doc],
        }

    def __str__(self):
        return '%s\t%s %s\t%s' % (self.id, self.type, self.arg, self.text)

    STANDOFF_RE = re.compile(r'^(\S+)\t(\S+) (\S+)\t(.*)$')

def parse_standoff_line(line):
    if not line:
        return None
    elif line[0] == 'T':
        return Textbound.from_standoff(line)
    elif line[0] == 'R':
        return Relation.from_standoff(line)
    elif line[0] == 'E':
        return Event.from_standoff(line)
    elif line[0] == 'N':
        return Normalization.from_standoff(line)
    elif line[0] in ('A', 'M'):
        return Attribute.from_standoff(line)
    elif line[0] == '#':
        return Comment.from_standoff(line)
    else:
        print >> sys.stderr, 'Warning: discarding unrecognized line:', line

def parse_standoff(source):
    annotations = []
    for line in source:
        line = line.rstrip('\n')
        if line.strip() == '':
            continue
        annotations.append(parse_standoff_line(line))
    return annotations

def get_source_db(files, options=None):
    if options is not None:
        return options.source
    else:
        # TODO: guess from filenames
        return DEFAULT_DB

def get_source_id(filenames, options=None):
    bases = [path.basename(fn) for fn in filenames]
    roots = [path.splitext(bn)[0] for bn in bases]
    uniques = set(roots)
    if len(uniques) > 1:
        print >> sys.stderr, 'Warning: ambiguous source: %s' % str(uniques)
    source = list(uniques)[0]
    # fix common variants
    m = re.match(r'^(?:(?:PMID|pubmed)[_-]?)?(\d+)$', source)
    if m:
        return m.group(1)
    else:
        return source

def to_pubannotation(annotations, text, files, options=None):
    ann_by_id = { a.id: a for a in annotations }
    pubann = defaultdict(list)
    for a in annotations:
        for category, annotations in a.to_pubannotation(ann_by_id).items():
            pubann[category].extend(annotations)
    pubann['sourcedb'] = get_source_db(files, options)
    pubann['sourceid'] = get_source_id(files, options)
    pubann['text'] = text
    return pubann

def pretty(doc):
    return json.dumps(doc, sort_keys=True, indent=2, separators=(',', ': '))

def find_texts(annotations, options=None):
    if options.textdir is not None:
        dirs = options.textdir
    else:
        dirs = [path.dirname(fn) for fn in annotations]
    base = path.basename(annotations[0])
    textbase = path.splitext(base)[0] + '.txt'
    candidates = [path.join(d, textbase) for d in dirs]
    return [fn for fn in candidates if path.exists(fn)]

def annotations_and_text(filenames, options=None):
    annotations, texts = [], []
    for fn in filenames:
        root, ext = path.splitext(path.basename(fn))
        if ext in ('.txt',):
            texts.append(fn)
        else:
            annotations.append(fn)
    if len(texts) == 0:
        texts = find_texts(annotations, options)
        assert texts != [], 'Failed to find text for %s' % str(annotations)
    if len(texts) > 1:
        print >> sys.stderr, 'Warning: multiple texts for %s' % str(annotations)
    return annotations, texts[0]

def verify_text(annotations, text):
    for a in annotations:
        a.verify_text(text)

def output_file_name(annotation, options):
    source = annotation.get('sourcedb')
    if source is None:
        prefix = ''
    elif source == 'PubMed':
        prefix = 'PMID-'
    else:
        print >> sys.stderr, 'Warning: unrecognized source %s' % source
        prefix = ''
    base = prefix + annotation['sourceid'] + '.json'
    return path.join(options.output, base)

def output_pubannotation(annotation, options=None):
    if options is None or options.output is None:
        print >> sys.stdout, pretty(annotation)
    else:
        with open(output_file_name(annotation, options), 'wt') as out:
            print >> out, pretty(annotation)

def process_files(files, options=None):
    annfiles, textfile = annotations_and_text(files, options)
    standoff = parse_standoff(fileinput.input(annfiles))
    text = codecs.open(textfile, encoding=DEFAULT_ENCODING).read()
    verify_text(standoff, text)
    pubann = to_pubannotation(standoff, text, files, options)
    output_pubannotation(pubann, options)

def group_files(filenames):
    files_by_basename = defaultdict(list)
    for fn in filenames:
        root, ext = path.splitext(path.basename(fn))
        files_by_basename[root].append(fn)
    return files_by_basename.values()

def main(argv):
    args = argparser().parse_args(argv[1:])

    for group in group_files(args.file):
        process_files(group, args)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
