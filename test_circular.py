import drawsvg as draw
from math import pi, cos, sin
from dataclasses import dataclass
from typing import List, Tuple
from mtSVG import Gene, MtGenome, Point, DrawableGenome, get_genomes, get_clean_name



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
RADIUS_RATIO = 0.88

RIBBON_HEIGHT = GENE_HEIGHT + ORIENTATION_HEIGHT + INTRA_GENOME_SPACE + INTER_GENOME_SPACE + SPECIES_HEIGHT

def get_color(color_scheme: dict, key: str) -> str:
    return next((color_scheme[k] for k in color_scheme if key.lower().startswith(k)), '#ffffff')


# UPDATED
def get_drawing(drawables: List[DrawableGenome], circular: False) -> draw.Drawing:
    width = max([(drawable.genome.get_scaled_length() * SCALE_FACTOR) + len(drawable.genome.genes) * STROKE_WIDTH
                 for drawable in drawables])
    if circular:
        width = int(width / pi)
        return draw.Drawing(width, width)
    height = RIBBON_HEIGHT * len(drawables)
    return draw.Drawing(width, height)

# NEW
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


# NEW
def x_to_deg(x: float, radius: float) -> float:
    angle = x / radius
    angle = pi/2 - angle if angle <= pi/2 else 5*pi/2 - angle
    return angle * (180 / pi)

# NEW
def x_to_polar(x: float, radius: float) -> Tuple[float, float]:
    angle = (x / radius) - pi/2
    return angle, radius * cos(angle), radius * sin(angle)
    
        
# NEW
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
    # TODO 
 
    return x_pos + gene.scaled_length * SCALE_FACTOR
    

# NEW
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


# ----------------------------- TESTING -----------------------------

gffs = [('Halocynthia roretzi', 14771, 'example/h_roretzi.gff', False)]
genomes = get_genomes(gffs, 'cox1', 0, False)
draw_circle(genomes, 'test.svg', monochromatic=False, font='Arial', full_name=False, oriented=False)
print('Done!')