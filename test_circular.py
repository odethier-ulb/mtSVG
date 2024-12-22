import drawsvg as draw
from math import pi
from dataclasses import dataclass
from typing import List, Tuple
from mtSVG import Gene, MtGenome, Point, DrawableGenome, get_genomes



# ----------------------------- DRAWING -----------------------------

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
PIXEL_SCALE = 0.3
INTER_GENOME_SPACE = 50
INTRA_GENOME_SPACE = 10
GENE_HEIGHT = 160
ORIENTATION_HEIGHT = 30
SPECIES_HEIGHT = 80

RIBBON_HEIGHT = GENE_HEIGHT + ORIENTATION_HEIGHT + INTRA_GENOME_SPACE + INTER_GENOME_SPACE + SPECIES_HEIGHT


# UPDATED
def get_drawing(drawables: List[DrawableGenome], circular: False) -> draw.Drawing:
    width = max([(drawable.genome.get_scaled_length() * SCALE_FACTOR) + len(drawable.genome.genes) * STROKE_WIDTH
                 for drawable in drawables])
    if circular:
        width = int(width / pi)
        return draw.Drawing(width, width)
    height = RIBBON_HEIGHT * len(drawables)
    return draw.Drawing(width, height)


def draw_circular_genome(drawable: DrawableGenome, drawing: draw.Drawing):
    # draw inner and outer circles
    c_x, c_y = drawing.width / 2, drawing.height / 2
    r_in, r_out = c_x - RIBBON_HEIGHT/2 - STROKE_WIDTH, c_x - STROKE_WIDTH
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
    

# NEW
def draw_circle(genomes: List[MtGenome], output: str,
                 monochromatic: bool = False,
                 font: str = 'Arial',
                 full_name: bool = False,
                 oriented: bool = False):
    
    color_scheme = [], COLOR_SCHEMES['monochromatic'] if monochromatic else COLOR_SCHEMES['default']
    drawable = DrawableGenome(Point(0, 0), color_scheme, font, full_name, oriented, genomes[0])
    drawing = get_drawing([drawable], circular=True)
    draw_circular_genome(drawable, drawing)
    drawing.set_pixel_scale(PIXEL_SCALE)
    drawing.save_svg(output)


# ----------------------------- TESTING -----------------------------

gffs = [('Styela plicata', 14414, 'example/s_plicata.gff', False)]
genomes = get_genomes(gffs, 'cox1', 0, False)
draw_circle(genomes, 'test.svg', monochromatic=False, font='Arial', full_name=False, oriented=False)
print('Done!')