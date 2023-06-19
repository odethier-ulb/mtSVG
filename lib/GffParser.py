import math
from typing import List, Tuple
from dataclasses import dataclass

GENE_CLASSES = set(['gene', 'trna', 'rrna'])


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


def __parse_gff(filepath: str) -> List[Gene]:
    genes = []
    with open(filepath, 'rt') as f:
        for line in f:
            lsplt = line.strip().split()
            if len(lsplt) < 9 or lsplt[2].lower() not in GENE_CLASSES:
                continue
            genes.append(Gene(lsplt[8].split('=')[-1].lower(), lsplt[6], int(lsplt[3]), int(lsplt[4])))
    return genes


def get_genomes(species: List[Tuple[str, int, str]], reverse: bool = False, start: str = 'cox1',
                intergenic: int = 100) -> List[MtGenome]:
    genomes, max_length = [MtGenome(sp[0], sp[1], __parse_gff(sp[2])) for sp in species], -1

    # reverse
    if reverse:
        for genome in genomes:
            genome.genes.reverse()

    # align to start gene
    for genome in genomes:
        start_idx = -1
        for i in range(0, len(genome.genes)):
            if start in genome.genes[i].name:
                start_idx = i
                break
        new_genes = []
        for i in range(start_idx, len(genome.genes)):
            new_genes.append(genome.genes[i])
            max_length = max(max_length, genome.genes[i].get_length(genome.length))
        for i in range(0, start_idx):
            new_genes.append(genome.genes[i])
            max_length = max(max_length, genome.genes[i].get_length(genome.length))
        genome.genes = new_genes

    # add intergenic regions
    if intergenic > 0:
        for genome in genomes:
            for i in range(len(genome.genes) - 1, 0, -1):
                if genome.genes[i].start < genome.genes[i-1].end < genome.genes[i].end:
                    # overlap
                    continue
                region_length = genome.genes[i].start - genome.genes[i-1].end
                if region_length < 0:
                    region_length = genome.genes[i].start + genome.length - genome.genes[i-1].end
                if region_length > intergenic:
                    genome.genes.insert(i, Gene('intergenic', None, genome.genes[i-1].end, genome.genes[i].start))
                    i -= 1

    # compute scaled length from min 1 to max 10
    unit = max_length / 10.
    for genome in genomes:
        for gene in genome.genes:
            gene.scaled_length = max(1, int(math.ceil(gene.get_length(genome.length) / unit)))

    return genomes
