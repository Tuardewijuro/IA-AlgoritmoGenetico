import random
from typing import List, Dict, Tuple

class Gene:
    def __init__(self, course: str, lecturer: str, room: str, day: str, time: str, prodi: str, sks: int, smt: int, student: int):
        self.course = course
        self.lecturer = lecturer
        self.room = room
        self.day = day
        self.time = time
        self.prodi = prodi
        self.sks = sks
        self.smt = smt
        self.student = student

class Chromosome:
    def __init__(self, genes: List[Gene]):
        self.genes = genes
        self.fitness = 0

class GeneticAlgorithm:
    def __init__(self, courses: List[Dict], rooms: List[str], days: List[str], times: List[str], 
                 population_size: int, mutation_rate: float, crossover_rate: float,
                 room_props: Dict[str, Dict], fitness_settings: Dict):
        self.courses = courses
        self.rooms = rooms
        self.days = days
        self.times = times
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.room_props = room_props
        self.fitness_settings = fitness_settings
        self.population = []

    def initialize_population(self):
        for _ in range(self.population_size):
            chromosome = Chromosome([
                Gene(course['course'], course['lecturer'], 
                     random.choice(self.rooms), 
                     random.choice(self.days), 
                     random.choice(self.times),
                     course['prodi'], course['sks'], course['smt'], course['student'])
                for course in self.courses
            ])
            self.population.append(chromosome)

    def fitness_function(self, chromosome: Chromosome) -> float:
        penalty = 0

        for gene in chromosome.genes:
            # Check room capacity
            room_capacity = self.room_props.get(gene.room, {}).get('capacity', 0)
            if room_capacity < gene.student:
                penalty += self.fitness_settings['roomOverCapacity']['penalty']
                
                # Check time constraints
                time_start = int(gene.time.split(' - ')[0].split('.')[0])
                if (time_start > 3 and time_start < 6) or (time_start < 4 and (time_start + gene.sks - 1) > 3):
                    penalty += 1
                
                if (gene.day == 'kamis' and time_start > 2 and time_start < 6) or \
                (gene.day == 'kamis' and time_start < 3 and (time_start + gene.sks - 1) > 2):
                    penalty += 1
        
        # Check for conflicts
        for i, gene1 in enumerate(chromosome.genes):
            for j, gene2 in enumerate(chromosome.genes[i+1:]):
                if gene1.day == gene2.day and gene1.time == gene2.time:
                    if gene1.room == gene2.room:
                        penalty += self.fitness_settings['sameRoomSameTime']['penalty']
                    if gene1.lecturer == gene2.lecturer:
                        penalty += self.fitness_settings['sameLecturerSameTime']['penalty']
                    if gene1.prodi == gene2.prodi and gene1.smt == gene2.smt:
                        penalty += self.fitness_settings['sameProdiSameSemesterSameTime']['penalty']
        
        return 1 / (1 + penalty)

    def tournament_selection(self, tournament_size: int) -> Chromosome:
        tournament = random.sample(self.population, tournament_size)
        return max(tournament, key=lambda chromo: chromo.fitness)

    def uniform_crossover(self, parent1: Chromosome, parent2: Chromosome) -> Tuple[Chromosome, Chromosome]:
        child1_genes = []
        child2_genes = []
        for gene1, gene2 in zip(parent1.genes, parent2.genes):
            if random.random() < 0.5:
                child1_genes.append(Gene(gene1.course, gene1.lecturer, gene1.room, gene1.day, gene1.time,
                                         gene1.prodi, gene1.sks, gene1.smt, gene1.student))
                child2_genes.append(Gene(gene2.course, gene2.lecturer, gene2.room, gene2.day, gene2.time,
                                         gene2.prodi, gene2.sks, gene2.smt, gene2.student))
            else:
                child1_genes.append(Gene(gene2.course, gene2.lecturer, gene2.room, gene2.day, gene2.time,
                                         gene2.prodi, gene2.sks, gene2.smt, gene2.student))
                child2_genes.append(Gene(gene1.course, gene1.lecturer, gene1.room, gene1.day, gene1.time,
                                         gene1.prodi, gene1.sks, gene1.smt, gene1.student))
        return Chromosome(child1_genes), Chromosome(child2_genes)

    def mutate(self, chromosome: Chromosome):
        for gene in chromosome.genes:
            if random.random() < self.mutation_rate:
                gene.room = random.choice(self.rooms)
                gene.day = random.choice(self.days)
                gene.time = random.choice(self.times)

    def evolve(self, generations: int):
        self.initialize_population()
        
        for _ in range(generations):
            new_population = []
            
            while len(new_population) < self.population_size:
                parent1 = self.tournament_selection(3)
                parent2 = self.tournament_selection(3)
                
                if random.random() < self.crossover_rate:
                    child1, child2 = self.uniform_crossover(parent1, parent2)
                else:
                    child1 = Chromosome([Gene(g.course, g.lecturer, g.room, g.day, g.time,
                                              g.prodi, g.sks, g.smt, g.student) for g in parent1.genes])
                    child2 = Chromosome([Gene(g.course, g.lecturer, g.room, g.day, g.time,
                                              g.prodi, g.sks, g.smt, g.student) for g in parent2.genes])
                
                self.mutate(child1)
                self.mutate(child2)
                
                for child in [child1, child2]:
                    child.fitness = self.fitness_function(child)
                
                new_population.extend([child1, child2])
            
            self.population = sorted(new_population, key=lambda x: x.fitness, reverse=True)[:self.population_size]
        
        best_solution = max(self.population, key=lambda chromo: chromo.fitness)
        return best_solution

# Input data
days = ["senin", "selasa", "rabu", "kamis", "jumat", "sabtu"]
times = ["08.00 - 08.50", "09.00 - 09.50", "10.00 - 10.50", "11.00 - 11.50", "19.00 - 19.50"]
rooms = ["R.1&2", "R.3", "R.4", "R.5", "R.6", "R.7", "R.8", "R.9&10"]
room_props = {
    "R.1&2": {"owner": ["ilmu hukum", "ilmu hukum (pagi)", "manajemen", "akuntansi"], "capacity": 150},
    "R.3": {"owner": ["all"], "capacity": 40},
    "R.4": {"owner": ["all"], "capacity": 40},
    "R.5": {"owner": ["all"], "capacity": 40},
    "R.6": {"owner": ["all"], "capacity": 40},
    "R.7": {"owner": ["all"], "capacity": 40},
    "R.8": {"owner": ["all"], "capacity": 40},
    "R.9&10": {"owner": ["all"], "capacity": 80}
}

room_props = {
    "R.1&2": {"owner": ["ilmu hukum", "ilmu hukum (pagi)", "manajemen", "akuntansi"], "capacity": 150}
}

courses = [
    {"prodi": "manajemen", "course": "pengantar ekonomi mikro (ma&ak)", "student": 30,
     "lecturer": "unna ria safitri", "sks": 3, "smt": 1, "required": True}
]

ga_params = {
    "nPopulations": 10000,
    "npops": 10,
    "nSelection": 3,
    "pCrossover": 0.85,
    "pMutation": 0.14,
    "fitnessThresshold": 0.8,
    "nSolution": 1,
    "stoppingCondition": {"nFitnessNoChange": 1000, "fitnessMax": 1}
}

fitness_settings = {
    "sameLecturerSameTime": {"enable": True, "penalty": 1},
    "sameProdiSameSemesterSameTime": {"enable": True, "penalty": 1},
    "sameRoomSameTime": {"enable": True, "penalty": 1},
    "timeOver": {"enable": True, "penalty": 1},
    "roomOverCapacity": {"enable": True, "penalty": 1},
    "roomUsedByOthers": {"enable": False, "penalty": 0.01},
    "sameLecturerSameDay": {"enable": False, "penalty": 0.001},
    "sameLecturerHasSequence": {"enable": False, "penalty": 0.001},
    "sameProdiSameSemesterSameDay": {"enable": False, "penalty": 0.001},
    "sameProdiSameSemesterHasSequence": {"enable": False, "penalty": 0.001}
}

# Initialize and run GA
ga = GeneticAlgorithm(
    courses=courses,
    rooms=rooms,
    days=days,
    times=times,
    population_size=ga_params["nPopulations"],
    mutation_rate=ga_params["pMutation"],
    crossover_rate=ga_params["pCrossover"],
    room_props=room_props,
    fitness_settings=fitness_settings
)

best_solution = ga.evolve(generations=ga_params["stoppingCondition"]["nFitnessNoChange"])

# Output
print({
    "i": 2632,
    "fitness": best_solution.fitness,
    "conflict": [],
    "schedule": [
        {
            "room": gene.room,
            "day": gene.day,
            "time": gene.time,
            "prodi": gene.prodi,
            "lecturer": gene.lecturer,
            "smt": gene.smt,
            "sks": gene.sks,
            "course": gene.course,
            "student": gene.student,
            "required": True
        } for gene in best_solution.genes
    ]
})