import drawsvg as draw
from typing import List
from dataclasses import dataclass
from lib.GffParser import MtGenome, Gene

COLOR_SCHEMES = {'default': {'co': '#f2ed8d', 'na': '#b6e07b', 'atp': '#b3e6e8',
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
PIXEL_SCALE = 1.5
INTER_GENOME_SPACE = 50
INTRA_GENOME_SPACE = 10
GENE_HEIGHT = 160
ORIENTATION_HEIGHT = 30
SPECIES_HEIGHT = 80

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
    genome: MtGenome


def get_color(color_scheme: dict, key: str) -> str:
    return next((color_scheme[k] for k in color_scheme if key.lower().startswith(k)), '#ffffff')


def get_drawables(genomes: List[MtGenome], color_scheme: dict, font: str) -> List[DrawableGenome]:
    drawables = []
    for i in range(len(genomes)):
        drawables.append(DrawableGenome(Point(0, i * RIBBON_HEIGHT), color_scheme, font, genomes[i]))
    return drawables


def get_drawing(drawables: List[DrawableGenome]) -> draw.Drawing:
    width = max([(drawable.genome.get_scaled_length() * SCALE_FACTOR) + len(drawable.genome.genes) * STROKE_WIDTH
                 for drawable in drawables])
    height = RIBBON_HEIGHT * len(drawables)
    return draw.Drawing(width, height)


def get_clean_name(gene_name: str) -> str:
    if gene_name is None or len(gene_name) < 1:
        return gene_name
    elif gene_name.lower().startswith('trn'):
        return 'trn' + gene_name[3].upper()
    else:
        return gene_name.split('_')[0].split('-')[0]


def draw_genome(drawable: DrawableGenome, drawing: draw.Drawing):
    # draw species
    species_font_size = SPECIES_HEIGHT * 0.75
    drawing.append(draw.Text(drawable.genome.species + f' ({drawable.genome.length:,} bp)', species_font_size,
                             drawable.origin.x,
                             drawable.origin.y + species_font_size,
                             font_family=drawable.font, font_style='italic', font_weight='bold'))
    # draw genes
    gene_origin = Point(drawable.origin.x + STROKE_WIDTH, drawable.origin.y + SPECIES_HEIGHT)
    for gene in drawable.genome.genes:
        gene_origin = draw_gene(gene, gene_origin, drawing, drawable.color_scheme, True, drawable.font)


def draw_gene(gene: Gene, origin: Point, drawing: draw.Drawing, color_scheme: dict, with_orient: bool, font: str) -> Point:
    # draw gene
    drawing.append(draw.Rectangle(origin.x, origin.y,
                                  gene.scaled_length * SCALE_FACTOR, GENE_HEIGHT,
                                  fill=get_color(color_scheme, gene.name), stroke='black', stroke_width=STROKE_WIDTH))

    # draw gene name
    gene_name = get_clean_name(gene.name)
    font_size = int(GENE_HEIGHT / 3)
    gene_size = len(gene_name) * font_size / 2

    if gene_size < gene.scaled_length * SCALE_FACTOR:
        drawing.append(draw.Text(gene_name, font_size,
                                 origin.x + ((gene.scaled_length * SCALE_FACTOR) - gene_size) / 2,
                                 origin.y + (GENE_HEIGHT + font_size * 0.75) / 2,
                                 font_family=font))
    else:
        if gene_size > GENE_HEIGHT:
            font_size = int(2 * GENE_HEIGHT / len(gene_name))
            gene_size = len(gene_name) * font_size / 2

        txt_x = origin.x + ((gene.scaled_length * SCALE_FACTOR) + font_size * 0.7) / 2
        txt_y = origin.y + (GENE_HEIGHT + gene_size) / 2
        drawing.append(draw.Text(gene_name, font_size, txt_x, txt_y,
                                 font_family=font, transform=f'rotate(270, {txt_x}, {txt_y})'))

    # draw orientation
    if with_orient and gene.name != 'intergenic':
        orientation_color = get_color(color_scheme, gene.orientation)
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
                                  close=False, fill=orientation_color))
    return Point(origin.x + gene.scaled_length * SCALE_FACTOR, origin.y)


def save_drawing(drawing: draw.Drawing, filepath: str):
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(filepath)


def draw_ribbons(genomes: List[MtGenome], filepath: str):
    drawables = get_drawables(genomes, COLOR_SCHEMES['default'], 'Arial')
    drawing = get_drawing(drawables)

    for drawable in drawables:
        draw_genome(drawable, drawing)

    save_drawing(drawing, filepath)
