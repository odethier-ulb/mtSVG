#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import argparse
import sys
import drawsvg as draw

gene_colors = {'co': '#f2ed8d', 'na' : '#b6e07b', 'rrn' : '#c7ace3', 'atp': '#b3e6e8', 'trn' : '#e69d97'}


def get_color(gene):
    for k in gene_colors:
        if gene.startswith(k):
            return gene_colors[k]
    return #ffffff


def compute_length(start, end, mtdna_size):
    length = end - start
    if length < 0:
        length = mtdna_size - start + end
    return length


def get_gene_list(gff_file, mtdna_size, start, reverse):
    result, max_length = [], -9999

    # get gene list and size
    with open(gff_file, 'rt') as f:
        for line in f:
            if not "\tgene" in line and not "\ttRNA" in line and not "\trRNA" in line:
                continue
            lsplt = line.strip().split()
            gene = lsplt[8].split('=')[-1]
            length = compute_length(int(lsplt[3]), int(lsplt[4]), mtdna_size)
            if length > max_length:
                max_length = length
            result.append([gene, length])

    if reverse == 'true':
        result.reverse()

    # change start and order
    start_idx = -1
    for i in range(0, len(result)):
        if start in result[i][0]:
            start_idx = i
            break
    new_result = []
    for i in range(start_idx, len(result)):
        new_result.append(result[i])
    for i in range(0, start_idx):
        new_result.append(result[i])
    result = new_result
    
    # compute graphical length
    unit = float(max_length/10)
    for entry in result:
        entry.append(max(1, int(math.ceil(entry[1]/unit))) * 100)
        
    return result


def add_rectangle(center, width, height, gene, drawing):
    rectangle = draw.Lines(center[0] - width/2, center[1] - height/2,
                           center[0] + width/2, center[1] - height/2,
                           center[0] + width/2, center[1] + height/2,
                           center[0] - width/2, center[1] + height/2,
                           close=True, fill=get_color(gene), stroke='black', stroke_width=10)
    drawing.append(rectangle)
    drawing.append(draw.Text(gene, 100, center[0] - 100, center[1] + 25))


def main(gff, svg, length, start, reverse):
    genes = get_gene_list(gff, length, start, reverse)

    if len(genes) == 0:
        sys.exit('Error : no gene found')
   
    width = sum([gene[2] for gene in genes])
    center = [-width/2 + genes[0][2] / 2, 0]

    d = draw.Drawing(width + 500, 500, origin='center')

    for i in range(len(genes)):
        gene = genes[i]
        add_rectangle(center, gene[2], 400, gene[0], d)
        if i < len(genes) - 1:
            center[0] += gene[2] / 2 + genes[i+1][2] / 2 

    d.set_pixel_scale(2)
    d.save_svg(svg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a MITOS2 mtDNA GFF to a linear SVG representation')
    parser.add_argument('--gff', type=str, help='The path to the gff file')
    parser.add_argument('--svg', type=str, help='The path of the SVG to create', default='linear_mtdna.svg')
    parser.add_argument('--length', type=int, help='The mtDNA length in bp (used to scale the output)')
    parser.add_argument('--start', type=str, help='Gene to use at the start of the ribbon', default='cox1')
    parser.add_argument('--reverse', type=str, help='If true, the gene order in the gff is reversed', default='false')
    args = parser.parse_args()

    if args.gff is None:
        sys.exit('Error : no gff file')

    if args.length is None:
        sys.exit('Error : no mtDNA length')

    from lib.GffParser import get_genomes

    genomes = get_genomes([('test', 14934, args.gff)])


    #main(args.gff, args.svg, args.length, args.start, args.reverse)
    #draw_test(args.svg)
    print('Done !')
