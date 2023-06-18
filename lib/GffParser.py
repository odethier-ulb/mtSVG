from dataclasses import dataclass


@dataclass
class Gene:
    name: str
    length: int
    scaled_length: int
    orientation: str


@dataclass
class MtGenome:
    species: str
    length: int
    genes: list

    def scaled_length(self, offset: int = 0) -> int:
        return sum([gene.scaled_length + 2 * offset for gene in self.genes])

