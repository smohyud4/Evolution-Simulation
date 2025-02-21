import os

# Parameters
populations = [25, 50, 75, 100]
mutations = [0, 0.01, 0.03, 0.05]
probabilities = [0, 0.1, 0.3, 0.5]
t_sizes = [2, 3, 4, 5]

def run_experiments(n):
    count = 0
    def create_command(pop, mut, prob, size, iteration):
        inter = r"C:\Users\Sulaiman Mohyuddin\anaconda3\python.exe"
        file_name = f"{pop}_{mut}_{prob}_{size}_{iteration}.csv"

        # Wrap the full command in escaped double quotes for Windows
        return f'cmd /c ""{inter}" lab2.py --n {pop} --p_m {mut} --p_c {prob} --trn_size {size} --csv_output "{file_name}""'


    for population in populations:
        for mutation in mutations:
            for probability in probabilities:
                for size in t_sizes:
                    for i in range(n):
                        command = create_command(population, mutation, probability, size, i+1)
                        os.system(command)  
                        count += 1
    print(count)

run_experiments(20)