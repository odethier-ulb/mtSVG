import drawsvg as draw
from typing import List
from dataclasses import dataclass
from lib.GffParser import MtGenome, Gene

COLOR_SCHEMES = {'default': {'co': '#f2ed8d', 'na': '#b6e07b', 'atp': '#b3e6e8',
                             'rrn': '#c7ace3',
                             'trn': '#e69d97',
                             'intergenic': '#000000'},
                 'monochromatic': {'co': '#ffffff', 'na': '#ffffff', 'atp': '#ffffff',
                                   'rrn': '#ffffff',
                                   'trn': '#ffffff',
                                   'intergenic': '#000000'}}

SCALE_FACTOR = 50
STROKE_WIDTH = 5
PIXEL_SCALE = 2
INTER_GENOME_SPACE = 50
INTRA_GENOME_SPACE = 10
GENE_HEIGHT = 100
ORIENTATION_HEIGHT = 20
SPECIES_HEIGHT = 30

RIBBON_HEIGHT = GENE_HEIGHT + ORIENTATION_HEIGHT + INTRA_GENOME_SPACE + INTER_GENOME_SPACE + SPECIES_HEIGHT


@dataclass
class Point:
    x: int
    y: int


@dataclass
class DrawableGenome:
    origin: Point
    color_scheme: dict
    genome: MtGenome


def get_gene_color(color_scheme, gene: Gene) -> str:
    return next((color_scheme[k] for k in color_scheme if gene.name.lower().startswith(k)), '#ffffff')


def get_drawables(genomes: List[MtGenome], color_scheme: dict) -> List[DrawableGenome]:
    drawables = []
    for i in range(len(genomes)):
        drawables.append(DrawableGenome(Point(0, i * RIBBON_HEIGHT), color_scheme, genomes[i]))
    return drawables


def get_drawing(drawables: List[DrawableGenome]) -> draw.Drawing:
    width = max([(drawable.genome[0].get_scaled_length() * SCALE_FACTOR) + len(drawable.genome[0].genes) * STROKE_WIDTH
                 for drawable in drawables])
    height = RIBBON_HEIGHT * len(drawables)
    return draw.Drawing(width, height)


def draw_genome(drawable: DrawableGenome, drawing: draw.Drawing):
    origin = drawable.origin
    for gene in drawable.genome[0].genes:
        origin = draw_gene(gene, origin, drawing, get_gene_color(drawable.color_scheme, gene))


def draw_gene(gene: Gene, origin: Point, drawing: draw.Drawing, color: str) -> Point:
    drawing.append(draw.Rectangle(origin.x, origin.y,
                                  gene.scaled_length * SCALE_FACTOR, GENE_HEIGHT,
                                  fill=color, stroke='black', stroke_width=STROKE_WIDTH))
    #     drawing = draw.Drawing(500, 500)
    #
    #     drawing.append(draw.Rectangle(10, 10, 300, 200, fill='#1248ff', stroke='black', stroke_width=10))
    #
    #     drawing.append(draw.Rectangle(10, 230, 240, 50, fill='#ffffff'))
    #
    #     drawing.append(draw.Lines(250, 230,
    #                               310, 255,
    #                               250, 280, close=False, fill='#ffffff'))
    #     drawing.save_svg(svg)
    return Point(origin.x + gene.scaled_length * SCALE_FACTOR + STROKE_WIDTH, origin.y)


def save_drawing(drawing: draw.Drawing, filepath: str):
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(filepath)


def draw_ribbons(genomes: List[MtGenome], filepath: str):
    drawables = get_drawables(genomes, COLOR_SCHEMES['default'])
    drawing = get_drawing(drawables)

    for drawable in drawables:
        draw_genome(drawable, drawing)

    save_drawing(drawing, filepath)
