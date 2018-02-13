#!/usr/bin/env python
'''Extract the reads mapping to a list of contigs from a bam file

All contigs not included in a given contig list are removed.

Fastq files are written to PREFIX.R1.fastq.gz, PREFIX.R2.fastq.gz, and
PREFIX.unpaired.fastq.gz.

This works similar to blobtools bamfilter but does not exclude unpaired reads.

Usage:
extract_csomes.py [-n NCHUNKS] [-f FORMAT] [-p THREADS] -o PREFIX -i TXTFILE
BAM

Options:
-o, --output_prefix=PREFIX    Output prefix
-f, --output_format=FORMAT    Output file format ('bam' or 'fastq')
[default:'bam']
-i, --include=TXTFILE         List of contig names to include (separated by
newlines)
-n, --nchunks=NCHUNKS         Number of chunks to split the region list. By
default, NCHUNKS=THREADS
-p, --threads=THREADS         Number of threads [default: 1]
'''
import logging as log
import os
try:
    from queue import Queue
except ImportError:
    from Queue import Queue
from threading import Thread
from uuid import uuid4

from docopt import docopt
from plumbum import (
    FG,
    local,
)
from plumbum.cmd import cat

from fmbiopy.fmbio import (
    merge_bams,
    to_fastq,
)

from fmbiopy.fmlist import split_into_chunks
from fmbiopy.fmlog import setup_log
from fmbiopy.fmsystem import capture_stdout

setup_log()


def start_workers(nchunks, queue, bam):
    '''Initialize the worker threads'''
    for i in range(nchunks):
        worker = BamExtractor(queue, bam)
        log.info('Starting BamExtractor thread %s', i)
        worker.daemon = True
        worker.start()


class BamExtractor(Thread):
    '''Worker thread for extracting contigs from bam'''

    def __init__(self, queue, bam):
        Thread.__init__(self)
        self.queue = queue
        self.bam = bam

    def run(self):
        '''Get the work from the queue and run samtools'''
        while True:
            contig_list, output_file = self.queue.get()
            contig_list = ' '.join(contig_list)
            log.info("Extracting contigs...")

            # For some reason I couldn't get plumbum to work here with threading
            command = ' '.join(['samtools view -bh', self.bam, contig_list,
                                '>', output_file])
            os.system(command)
            log.info("Done extracting contigs, shutting down...")
            self.queue.task_done()


def main(bam, contig_file, output_format, nthreads, output_prefix, nchunks):
    log.info('Running extract_csomes.py with parameters:')
    log.info('Input Bam: %s', str(bam))
    log.info('Contig File: %s', str(contig_file))
    log.info('Output Format: %s', str(output_format))
    log.info('Threads: %s', str(nthreads))
    log.info('Output Prefix: %s', str(output_prefix))
    log.info('Number of chunks: %s', str(nchunks))
    with local.tempdir() as tmpdir:

        log.info('Reading input contig list')
        region_list = capture_stdout(cat[contig_file])
        if len(region_list) < nchunks:
            nchunks = len(region_list)

        log.info('Dividing contig list into %s chunks', nchunks)
        chunks = split_into_chunks(region_list, nchunks)
        chunk_bams = [tmpdir / (uuid4().hex + '.bam') for _ in range(nchunks)]

        queue = Queue(maxsize=nthreads)
        start_workers(nchunks, queue, bam)

        log.info('Queueing jobs')
        for chunk, chunk_bam in zip(chunks, chunk_bams):
            queue.put((chunk, chunk_bam))

        queue.join()

        if output_format == 'bam':
            output_bam = local.path(output_prefix + '.bam')
        else:
            output_bam = tmpdir / 'merged.bam'

        log.info('Merging temporary .bam files to %s', output_bam)
        merge_bams(chunk_bams, output_bam, sort_by="name")

        if output_format == 'fastq':
            log.info('Converting to .fastq')
            to_fastq(output_bam, output_prefix)


if __name__ == '__main__':
    opts = docopt(__doc__)
    nthreads = int(opts['--threads'])
    if not opts['--nchunks']:
        nchunks = nthreads
    else:
        nchunks = int(opts['--nchunks'])
    if nchunks < nthreads:
        raise ValueError("Number of chunks must be >= nthreads")

    output_format = opts['--output_format']
    if output_format not in ['bam', 'fastq']:
        raise ValueError("Unrecognized output format")

    main(bam=local.path(opts['BAM']),
         contig_file=local.path(opts['--include']),
         output_format=output_format,
         nthreads=nthreads,
         output_prefix=opts['--output_prefix'],
         nchunks=nchunks)
