#!/usr/bin/env python

import sys
import re
import fileinput
import json

from collections import defaultdict
from itertools import count
from os import path

def argparser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', metavar='DIR', default=None,
                        help='Output directory')
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

    STANDOFF_RE = None

    def to_pubannotation(self, ann_by_id):
        raise NotImplementedError

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

    def get_spans(self):
        spans = []
        for span in self.spans.split(';'):
            start, end = span.split(' ')
            spans.append((int(start), int(end)))
        return spans

    def to_pubannotation(self, ann_by_id):        
        spans = self.get_spans()
        start, end = spans[0][0], spans[-1][1]
        if len(spans) > 1:
            print >> sys.stderr, 'Warning: flattening span %s to %d-%d' % \
                (self.spans, start, end)
        doc = {
            'id': self.id,
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

    def get_args(self):
        a1, a2 = self.args.split(' ')
        a1key, a1val = a1.split(':', 1)
        a2key, a2val = a2.split(':', 1)
        return ((a1key, a1val), (a2key, a2val))
    
    def to_pubannotation(self, ann_by_id):
        arg1, arg2 = self.get_args()
        doc = {
            'id': self.id,
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

    def get_args(self):
        return [a.split(':', 1) for a in self.args.split(' ')]        

    def to_pubannotation(self, ann_by_id):
        relations = []
        for key, val in self.get_args():
            rid = new_id('R', ann_by_id)
            ann_by_id[rid] = None # reserve
            doc = {
                'id': rid,
                'pred': 'has'+key,
                'subj': val,
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

    def to_pubannotation(self, ann_by_id):
        # map to denotation + relation using Pierre Zweigenbaum's approach:
        # https://github.com/linkedannotation/blah2015/wiki/PubAnnotation-system
        spans = ann_by_id[self.arg].get_spans()
        start, end = spans[0][0], spans[-1][1]
        denotation = {
            'id': self.id,
            'obj': self.ref,
            'span': { 'begin': start, 'end': end },
        }
        rid = new_id('R', ann_by_id)
        ann_by_id[rid] = None # reserve
        relation = {
            'id': rid,
            'pred': 'Normalization',
            'subj': self.arg,
            'obj': self.id,
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

    def to_pubannotation(self, ann_by_id):
        pred = self.type + (self.val if self.val is not None else '')
        doc = {
            'id': self.id,
            'pred': pred,
            'obj': self.arg,
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

    def to_pubannotation(self, ann_by_id):
        print >> sys.stderr, 'Warning: comment conversion TODO'
        return {}

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

def convert_standoff(annotations):
    ann_by_id = { a.id: a for a in annotations }
    pubann = defaultdict(list)
    for a in annotations:
        for category, annotations in a.to_pubannotation(ann_by_id).items():
            pubann[category].extend(annotations)
    print pubann
    print pretty(pubann)
    return pubann

def pretty(doc):
    return json.dumps(doc, sort_keys=True, indent=2, separators=(',', ': '))

def process_files(files, options=None):
    standoff = parse_standoff(fileinput.input(files))
    pubann = convert_standoff(standoff)
    print pretty(pubann)

def group_files(filenames):
    files_by_basename = defaultdict(list)
    for fn in filenames:
        bn = path.basename(fn)
        root, ext = path.splitext(bn)
        if ext in ('.txt',):
            print >> sys.stderr, 'Note: not processing text file', fn
        else:
            files_by_basename[root].append(fn)
    return files_by_basename.values()

def main(argv):
    args = argparser().parse_args(argv[1:])

    for group in group_files(args.file):
        process_files(group)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
