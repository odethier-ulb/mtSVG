#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse
import sys
import logging
import drawsvg as draw
from math import ceil, pi, cos, sin
from typing import List, Tuple
from dataclasses import dataclass


logging.basicConfig(
    level=logging.WARNING,  
    format="%(asctime)s - %(levelname)s - %(message)s",  
)


# ----------------------------- GFF PARSING -----------------------------

GENE_CLASSES_MITOS = set(['gene', 'trna', 'rrna'])
GENE_CLASSES_GENEBANK = set(['gene'])

@dataclass
class Gene:
    name: str
    orientation: str
    start: int
    end: int
    scaled_length: int = None

    def get_length(self, genome_size) -> int:
        length = self.end - self.start
        if length < 0:
            length = genome_size - self.start + self.end
        return length


@dataclass
class MtGenome:
    species: str
    length: int
    genes: list

    def get_scaled_length(self) -> int:
        return sum([g.scaled_length for g in self.genes])


def product_to_gene_name(product: str) -> str:
    p = product.lower()
    if '16s' in p or 'large' in p:
        return 'rrnl'
    elif '12s' in p or 'small' in p:
        return 'rrns'
    elif 'trna-ala' in p:
        return 'trnA'
    elif 'trna-arg' in p:
        return 'trnR'
    elif 'trna-asn' in p:
        return 'trnN'
    elif 'trna-asp' in p:
        return 'trnD'
    elif 'trna-asx' in p:
        return 'trnB'
    elif 'trna-cys' in p:
        return 'trnC'
    elif 'trna-gln' in p:
        return 'trnQ'
    elif 'trna-glu' in p:
        return 'trnE'
    elif 'trna-glx' in p:
        return 'trnZ'
    elif 'trna-gly' in p:
        return 'trnG'
    elif 'trna-his' in p:
        return 'trnH'
    elif 'trna-ile' in p:
        return 'trnI'
    elif 'trna-leu' in p:
        return 'trnL'
    elif 'trna-lys' in p:
        return 'trnK'
    elif 'trna-met' in p:
        return 'trnM'
    elif 'trna-phe' in p:
        return 'trnF'
    elif 'trna-pro' in p:
        return 'trnP'
    elif 'trna-ser' in p:
        return 'trnS'
    elif 'trna-thr' in p:
        return 'trnT'
    elif 'trna-trp' in p:
        return 'trnW'
    elif 'trna-tyr' in p:
        return 'trnY'
    elif 'trna-val' in p:
        return 'trnV'
    else: 
        raise Exception(f'Unknown product name: {p}')


def parse_gff(filepath: str, to_skip: List[str]) -> List[Gene]:
    genes = []
    with open(filepath, 'rt') as f:
        lines = f.readlines()
    lines = [line for line in lines if not line.startswith('#') and len(line.split()) >= 9]
    for idx, line in enumerate(lines):
        lsplt = line.strip().split('\t')
        is_mitos = lsplt[1].lower().startswith('mit')
        if is_mitos and lsplt[2].lower() in GENE_CLASSES_MITOS:
            gene_name = next((x.strip() for x in lsplt[8].split(";") if x.strip().startswith("Name=")), None)
            genes.append(Gene(gene_name.split('=')[-1].lower(), lsplt[6], int(lsplt[3]), int(lsplt[4])))
        elif not is_mitos and lsplt[2].lower() in GENE_CLASSES_GENEBANK:
            gene_name = next((x.strip() for x in lsplt[8].split(";") if x.strip().startswith("gene=")), None)
            if not gene_name is None:
                gene_name = gene_name.split('=')[-1].lower()
            else: 
                gene_name = next((x.strip() for x in lines[idx+1].strip().split()[8].split(";") if x.strip().startswith("product=")), None) 
                gene_name = product_to_gene_name(gene_name.split('=')[-1])
            genes.append(Gene(gene_name, lsplt[6], int(lsplt[3]), int(lsplt[4])))
        else:
            continue
    # filter out genes to skip and order by start
    genes_filtered = []
    for gene in genes:
        if not any([gene.name.startswith(s) for s in to_skip]):
            genes_filtered.append(gene)
    genes_filtered.sort(key=lambda x: x.start)
    return genes_filtered


def get_genomes(species: List[Tuple[str, int, str, bool]], start: str, intergenic: int, linear: bool, to_skip: str) -> List[MtGenome]:
    to_skip = [] if to_skip is None else [s.strip() for s in to_skip.split(',')]
    tmp_genomes, genomes, max_length = [MtGenome(sp[0], sp[1], parse_gff(sp[2], to_skip)) for sp in species], [], -1

    # filter out genomes with no genes
    for genome in tmp_genomes:
        if len(genome.genes) > 0:
            genomes.append(genome)
        else:
            logging.warning(f'No gene found for {genome.species}, removed from the drawing')

    # stop if no genome
    if len(genomes) == 0:
        sys.exit('Error : no gene found in any genome')

    # add intergenic regions
    if intergenic > 0:
        for genome in genomes:
            for i in range(len(genome.genes) - 1, 0, -1):
                if genome.genes[i].start < genome.genes[i - 1].end < genome.genes[i].end:
                    # overlap
                    continue
                region_length = genome.genes[i].start - genome.genes[i - 1].end
                if region_length < 0:
                    region_length = genome.genes[i].start + genome.length - genome.genes[i - 1].end
                if region_length >= intergenic:
                    genome.genes.insert(i, Gene('intergenic', None, genome.genes[i - 1].end, genome.genes[i].start))
                    i -= 1
            # fill the gap between last gene and total length
            if genome.genes[-1].end < genome.length:
                genome.genes.append(Gene('intergenic', None, genome.genes[-1].end, genome.length))
    # reverse
    for i in range(len(species)):
        if species[i][3]:
            genomes[i].genes.reverse()

    # align to start gene
    for genome in genomes:
        start_idx = -1
        if linear and genome.length > 0:
            start = genome.genes[0].name
        for i in range(0, len(genome.genes)):
            if start.lower() in genome.genes[i].name.lower():
                start_idx = i
                break
        if start_idx == -1:
            logging.warning(f'Start gene {start} not found in {genome.species}, will use first gene found')
            start_idx = 0
        new_genes = []
        for i in range(start_idx, len(genome.genes)):
            new_genes.append(genome.genes[i])
            max_length = max(max_length, genome.genes[i].get_length(genome.length))
        for i in range(0, start_idx):
            new_genes.append(genome.genes[i])
            max_length = max(max_length, genome.genes[i].get_length(genome.length))
        genome.genes = new_genes

    # compute scaled length from min 1 to max 10
    unit = max_length / 10.
    for genome in genomes:
        for gene in genome.genes:
            gene.scaled_length = max(1, int(ceil(gene.get_length(genome.length) / unit)))

    return genomes


# ----------------------------- DRAWING -----------------------------

COLOR_SCHEMES = {'default': {'co': '#f2ed8d', 'cy': '#f2ed8d',
                             'na': '#b6e07b', 'nd': '#b6e07b', 
                             'atp': '#b3e6e8',
                             'rrn': '#c7ace3',
                             'trn': '#e69d97',
                             'intergenic': '#000000',
                             '+': '#a8d2e7', '-': '#ac759a'},
                 'monochromatic': {'co': '#ffffff', 'na': '#ffffff', 'atp': '#ffffff',
                                   'rrn': '#ffffff',
                                   'trn': '#ffffff',
                                   'intergenic': '#000000',
                                   '+': '#000000', '-': '#000000'}}

SCALE_FACTOR = 50
STROKE_WIDTH = 5
PIXEL_SCALE = 0.3
INTER_GENOME_SPACE = 50
INTRA_GENOME_SPACE = 10
GENE_HEIGHT = 160
ORIENTATION_HEIGHT = 30
SPECIES_HEIGHT = 80
RADIUS_RATIO = 0.88

RIBBON_HEIGHT = GENE_HEIGHT + ORIENTATION_HEIGHT + INTRA_GENOME_SPACE + INTER_GENOME_SPACE + SPECIES_HEIGHT


@dataclass
class Point:
    x: int
    y: int


@dataclass
class DrawableGenome:
    origin: Point
    color_scheme: dict
    font: str
    full_name: bool
    oriented: bool
    genome: MtGenome


def get_color(color_scheme: dict, key: str) -> str:
    return next((color_scheme[k] for k in color_scheme if key.lower().startswith(k)), '#ffffff')


def get_drawing(drawables: List[DrawableGenome], circular=False) -> draw.Drawing:
    width = max([(drawable.genome.get_scaled_length() * SCALE_FACTOR) + len(drawable.genome.genes) * STROKE_WIDTH
                 for drawable in drawables])
    if circular:
        width = int(width / pi)
        return draw.Drawing(width, width)
    height = RIBBON_HEIGHT * len(drawables)
    return draw.Drawing(width, height)


def get_clean_name(gene_name: str) -> str:
    try:
        if gene_name.lower().startswith('trn'):
            return 'trn' + gene_name[3].upper()
        else:
            return gene_name.split('_')[0].split('-')[0]
    except:
        return gene_name


#----------------------------- RIBBON -----------------------------#

def draw_genome(drawable: DrawableGenome, drawing: draw.Drawing):
    # draw species
    species_font_size = SPECIES_HEIGHT * 0.75
    drawing.append(draw.Text(drawable.genome.species + f' ({drawable.genome.length:,} bp)', species_font_size,
                             drawable.origin.x,
                             drawable.origin.y + species_font_size,
                             font_family=drawable.font,
                             font_style='italic',
                             font_weight='bold'))
    # draw genes
    gene_origin = Point(drawable.origin.x + STROKE_WIDTH, drawable.origin.y + SPECIES_HEIGHT)
    for gene in drawable.genome.genes:
        gene_origin = draw_gene(drawable, gene, gene_origin, drawing)


def draw_gene(drawable: DrawableGenome, gene: Gene, origin: Point, drawing: draw.Drawing) -> Point:
    # draw gene
    drawing.append(draw.Rectangle(origin.x, origin.y,
                                  gene.scaled_length * SCALE_FACTOR, GENE_HEIGHT,
                                  fill=get_color(drawable.color_scheme, gene.name),
                                  stroke='black',
                                  stroke_width=STROKE_WIDTH))

    # draw gene name
    gene_name = gene.name if drawable.full_name else get_clean_name(gene.name)
    font_size = int(GENE_HEIGHT / 3)
    gene_size = len(gene_name) * font_size / 2

    if gene_size < gene.scaled_length * SCALE_FACTOR:
        drawing.append(draw.Text(gene_name, font_size,
                                 origin.x + ((gene.scaled_length * SCALE_FACTOR) - gene_size) / 2,
                                 origin.y + (GENE_HEIGHT + font_size * 0.75) / 2,
                                 font_family=drawable.font))
    else:
        if gene_size > GENE_HEIGHT:
            font_size = int(2 * GENE_HEIGHT / len(gene_name))
            gene_size = len(gene_name) * font_size / 2

        txt_x = origin.x + ((gene.scaled_length * SCALE_FACTOR) + font_size * 0.7) / 2
        txt_y = origin.y + (GENE_HEIGHT + gene_size) / 2
        drawing.append(draw.Text(gene_name, font_size, txt_x, txt_y,
                                 font_family=drawable.font,
                                 transform=f'rotate(270, {txt_x}, {txt_y})'))

    # draw orientation
    if drawable.oriented and gene.name != 'intergenic':
        orientation_color = get_color(drawable.color_scheme, gene.orientation)
        origin_x = origin.x if gene.orientation == '+' else origin.x + SCALE_FACTOR
        if gene.scaled_length > 1:
            drawing.append(draw.Rectangle(origin_x, origin.y + GENE_HEIGHT + INTRA_GENOME_SPACE,
                                          (gene.scaled_length - 1) * SCALE_FACTOR, ORIENTATION_HEIGHT,
                                          fill=orientation_color))
            origin_x += (gene.scaled_length - 1) * SCALE_FACTOR if gene.orientation == '+' else 0
        # draw arrow
        arrow_x = origin_x + SCALE_FACTOR if gene.orientation == '+' else origin_x - SCALE_FACTOR
        drawing.append(draw.Lines(origin_x, origin.y + GENE_HEIGHT + INTRA_GENOME_SPACE,
                                  arrow_x, origin.y + GENE_HEIGHT + INTRA_GENOME_SPACE + ORIENTATION_HEIGHT / 2,
                                  origin_x, origin.y + GENE_HEIGHT + INTRA_GENOME_SPACE + ORIENTATION_HEIGHT,
                                  close=False,
                                  fill=orientation_color))
    return Point(origin.x + gene.scaled_length * SCALE_FACTOR, origin.y)


def draw_ribbons(genomes: List[MtGenome], output: str,
                 monochromatic: bool = False,
                 font: str = 'Arial',
                 full_name: bool = False,
                 oriented: bool = False):
    drawables, color_scheme = [], COLOR_SCHEMES['monochromatic'] if monochromatic else COLOR_SCHEMES['default']
    for i in range(len(genomes)):
        drawables.append(DrawableGenome(Point(0, i * RIBBON_HEIGHT),
                                        color_scheme, font, full_name, oriented,
                                        genomes[i]))
    drawing = get_drawing(drawables)
    for drawable in drawables:
        draw_genome(drawable, drawing)
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(output)


# ----------------------------- CIRCLE -----------------------------#

def x_to_deg(x: float, radius: float) -> float:
    angle = x / radius
    angle = pi/2 - angle if angle <= pi/2 else 5*pi/2 - angle
    return angle * (180 / pi)


def x_to_polar(x: float, radius: float) -> Tuple[float, float]:
    angle = (x / radius) - pi/2
    return angle, radius * cos(angle), radius * sin(angle)


def draw_circular_genome(drawable: DrawableGenome, drawing: draw.Drawing):
    # draw inner and outer circles
    c_x, c_y = drawing.width / 2, drawing.height / 2
    r_out = (drawable.genome.get_scaled_length() * SCALE_FACTOR) / (pi * 2)
    r_in = r_out - RIBBON_HEIGHT / 2
    drawing.append(draw.Circle(c_x, c_y, r_in, fill='none', stroke_width=STROKE_WIDTH, stroke='black'))
    drawing.append(draw.Circle(c_x, c_y, r_out, fill='none', stroke_width=STROKE_WIDTH, stroke='black'))
    # draw species name
    species_font_size = SPECIES_HEIGHT * 0.75
    sp_name, sp_length = drawable.genome.species, f'({drawable.genome.length:,} bp)'
    drawing.append(draw.Text(sp_name, species_font_size,
                             c_x - (len(sp_name) * species_font_size / 2) / 2, c_y - (species_font_size * 0.7) / 2,
                             font_family=drawable.font, font_style='italic', font_weight='bold'))
    drawing.append(draw.Text(sp_length, species_font_size,
                             c_x - (len(sp_length) * species_font_size / 2) / 2, c_y + (species_font_size * 0.7),
                             font_family=drawable.font, font_style='italic', font_weight='bold'))
    # draw genes
    x_pos = 0
    for gene in drawable.genome.genes:
        x_pos = draw_circular_gene(drawable, gene, c_x, c_y, x_pos, r_in, r_out, drawing)

        
def draw_circular_gene(drawable: DrawableGenome, gene: Gene, c_x: float, c_y: float, 
                       x_pos: float, r_in: float, r_out: float, drawing: draw.Drawing) -> float:
    # draw gene arcs
    angle_from = x_to_deg(x_pos, r_out)
    angle_to = x_to_deg(x_pos + gene.scaled_length * SCALE_FACTOR, r_out)
    drawing.append(draw.ArcLine(c_x, c_y, (r_in + r_out) / 2, angle_to, angle_from,
        stroke='black', stroke_width=RIBBON_HEIGHT/2 - STROKE_WIDTH, fill='none', fill_opacity=0.0))
    drawing.append(draw.ArcLine(c_x, c_y, (r_in + r_out) / 2, angle_to + 0.15, angle_from - 0.15,
        stroke=get_color(drawable.color_scheme, gene.name), stroke_width=RIBBON_HEIGHT/2 - STROKE_WIDTH, fill='none', fill_opacity=0.0))

    # draw gene name
    gene_name = gene.name if drawable.full_name else get_clean_name(gene.name)
    font_size = int(GENE_HEIGHT / 3.5)
    gene_size = len(gene_name) * font_size / 2

    if gene_size < gene.scaled_length * SCALE_FACTOR:
        angle, origin_x, origin_y = x_to_polar(x_pos + ((gene.scaled_length * SCALE_FACTOR) - gene_size) / 2, r_out)
        text_rotation = int((angle + pi/1.9) * (180/pi))
        text_x, text_y = RADIUS_RATIO * origin_x + c_x, RADIUS_RATIO * origin_y + c_y
        drawing.append(draw.Text(gene_name, font_size,text_x, text_y, font_family=drawable.font, 
                                 transform=f'rotate({text_rotation}, {text_x}, {text_y})'))   
    else:
        if gene_size > GENE_HEIGHT:
            font_size = int(2 * GENE_HEIGHT / len(gene_name))
            gene_size = len(gene_name) * font_size / 2

        angle, origin_x, origin_y = x_to_polar(x_pos + ((gene.scaled_length * SCALE_FACTOR) + font_size * 0.7) / 2, r_out)
        text_x, text_y = (RADIUS_RATIO - .03) * origin_x + c_x, (RADIUS_RATIO - .03) * origin_y + c_y
        text_rotation = int(angle * (180/pi))
        drawing.append(draw.Text(gene_name, font_size,text_x, text_y, font_family=drawable.font, 
                                 transform=f'rotate({text_rotation}, {text_x}, {text_y})'))

    # draw orientation
    if drawable.oriented and gene.name != 'intergenic':
        orientation_color = get_color(drawable.color_scheme, gene.orientation)
        origin_x = x_pos if gene.orientation == '+' else x_pos + SCALE_FACTOR
        angle_from = x_to_deg(origin_x, r_out)
        angle_to = x_to_deg(x_pos + (gene.scaled_length - 1) * SCALE_FACTOR, r_out)
        r_orientation = r_out - RIBBON_HEIGHT/1.8

        if gene.scaled_length > 1:
            drawing.append(draw.ArcLine(c_x, c_y, r_orientation, angle_to, angle_from, 
                                        stroke=orientation_color, stroke_width=ORIENTATION_HEIGHT / 2, 
                                        fill='none', fill_opacity=0.0))
            origin_x += (gene.scaled_length - 1) * SCALE_FACTOR if gene.orientation == '+' else 0

        # draw arrow
        arrow_x = origin_x + SCALE_FACTOR if gene.orientation == '+' else origin_x - SCALE_FACTOR
        _, x_1, y_1 = x_to_polar(origin_x, r_out)
        _, x_2, y_2 = x_to_polar(arrow_x, r_out)
        _, x_3, y_3 = x_to_polar(origin_x, r_out)
        ratio_1 = (r_orientation + ORIENTATION_HEIGHT / 2) / r_out
        ratio_2 = r_orientation / r_out
        ratio_3 = (r_orientation - ORIENTATION_HEIGHT / 2) / r_out
        drawing.append(draw.Lines(x_1 * ratio_1 + c_x, y_1 * ratio_1 + c_y, 
                                  x_2 * ratio_2 + c_x, y_2 * ratio_2  + c_y, 
                                  x_3 * ratio_3 + c_x, y_3 * ratio_3 + c_y, 
                                  close=False, fill=orientation_color))
 
    return x_pos + gene.scaled_length * SCALE_FACTOR
    

def draw_circle(genomes: List[MtGenome], output: str,
                 monochromatic: bool = False,
                 font: str = 'Arial',
                 full_name: bool = False, 
                 oriented: bool = False):
    
    color_scheme = COLOR_SCHEMES['monochromatic'] if monochromatic else COLOR_SCHEMES['default']
    drawable = DrawableGenome(Point(0, 0), color_scheme, font, full_name, oriented, genomes[0])
    drawing = get_drawing([drawable], circular=True)
    draw_circular_genome(drawable, drawing)
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(output)


# ----------------------------- MAIN -----------------------------#

def parse_gffs(filepath: str) -> List[Tuple[str, int, str, bool]]:
    try:
        results = []
        with open(filepath, 'rt') as f:
            for line in f:
                lstrip = line.strip()
                if lstrip.startswith('#'):
                    continue
                lsplt = lstrip.split(';')
                if len(lsplt) < 4:
                    continue
                results.append((lsplt[0], int(lsplt[1]), lsplt[2], bool(lsplt[3])))
        return results
    except:
        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a mtDNA GFF to a linear SVG representation')
    parser.add_argument('--gff', type=str, help='The path of a single gff file')
    parser.add_argument('--species', type=str, help='The species name (ignored if --gffs is used)')
    parser.add_argument('--size', type=int, help='The size of the mtDNA in base pair (ignored if --gffs is used)')
    parser.add_argument('--reversed', action='store_true', help='Reverse the gene order (ignored if --gffs is used')
    parser.add_argument('--gffs', type=str,
                        help='The path of the semicolon separated config file to draw multiple ribbons. Each entry in '
                             'the config file must have the following format: "species;mtdna size;gff path;reverse '
                             'gene order" like for instance "Halocynthia roretzi;14771;example/h_roretzi.gff;false". '
                             'Comments can be inserted using "#" to start a line and the last "true/false" value for '
                             'the gene order can be omitted when using false.')
    parser.add_argument('--start', type=str, help='Start gene of the ribbon', default='cox1')
    parser.add_argument('--linear', action='store_true',
                        help='Show the genes in the same order as in the gff file.')
    parser.add_argument('--intergenic', type=int,
                        help='Display intergenic regions having a size = or >, skipped if 0 is set', default=0)
    parser.add_argument('--oriented', action='store_true', help='Display gene orientations')
    parser.add_argument('--full_name', action='store_true', help='Display gene full names')
    parser.add_argument('--monochromatic', action='store_true', help='Do not colorize')
    parser.add_argument('--circular', action='store_true', help='Draw a circular representation (for single gff only)')
    parser.add_argument('--font', type=str, help='The font to use', default='Arial')
    parser.add_argument('--output', type=str, help='The path of the output to create', default='mtDNA.svg')
    parser.add_argument('--skip', type=str, help='Comma separated list of gene names to skip')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.gff is not None:
        if args.species is None:
            sys.exit('Error : missing species')
        if args.size is None:
            sys.exit('Error : missing size')
        gffs = [(args.species, args.size, args.gff, args.reversed)]
    elif args.gffs is not None:
        gffs = parse_gffs(args.gffs)
        if gffs is None:
            sys.exit('Error : wrong gffs file format')
        elif args.circular:
            sys.exit('Error : circular representation not supported with multiple genomes')
    else:
        sys.exit('Error : missing gff(s) file')

    genomes = get_genomes(gffs, args.start, args.intergenic, args.linear, args.skip)
    if args.circular:
        draw_circle(genomes, args.output, args.monochromatic, args.font, args.full_name, args.oriented)
    else:
        draw_ribbons(genomes, args.output, args.monochromatic, args.font, args.full_name, args.oriented)
    print('Done !')
