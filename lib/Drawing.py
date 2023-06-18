import drawsvg as draw
from dataclasses import dataclass
from GffParser import MtGenome, Gene

COLOR_SCHEMES = {'default': {'co': '#f2ed8d', 'na': '#b6e07b', 'rrn': '#c7ace3', 'atp': '#b3e6e8', 'trn': '#e69d97'}}
STROKE_WIDTH = 5
PIXEL_SCALE = 2
INTER_GENOME_SPACE = 50
INTRA_GENOME_SPACE = 10
GENE_HEIGHT = 100
ORIENTATION_HEIGHT = 20
SPECIES_HEIGHT = 30


@dataclass
class Point:
    x: int
    y: int


@dataclass
class DrawableGenome:
    origin: Point
    color_scheme: dict
    genome: MtGenome


def get_drawing(drawables: list[DrawableGenome]) -> draw.Drawing:
    width = max([drawable.genome.scaled_length(STROKE_WIDTH) for drawable in drawables])
    height = (GENE_HEIGHT + ORIENTATION_HEIGHT + INTRA_GENOME_SPACE + INTER_GENOME_SPACE + SPECIES_HEIGHT) * len(drawables)
    return draw.Drawing(width, height)


def draw_gene(gene: Gene, origin: Point, drawing: draw.Drawing) -> Point:
    return None


def save_drawing(drawing: draw.Drawing, filepath: str):
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(filepath)
