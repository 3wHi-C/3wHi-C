import random
file_path = "/mnt/disk1/1/gy/HiPore-C/distance/HiPoreC-K562-all.txt"
data_dict = {}
index_counts = {}
with open(file_path, "r") as file:
    for line in file:
        columns = line.strip().split("\t")
        index = columns[0]
        if index not in data_dict:
            data_dict[index] = []
            index_counts[index] = 0
        data_dict[index].append(line.strip())
        index_counts[index] += 1

for index, count in index_counts.items():
    if count >= 3:
        file_names = ["HiPoreC_random_0.txt","HiPoreC_random_1.txt","HiPoreC_random_2.txt"]
        selected_data = data_dict[index]
        #random.shuffle(selected_data)
        num_of_iterations = count-2
        for i in range(num_of_iterations):
            selected_data_1 = selected_data[i]
            selected_data_2 = selected_data[i+1]
            selected_data_3 = selected_data[i+2]
            with open(file_names[0],"a") as file:
                file.write(selected_data_1 + "\n")    
            with open(file_names[1],"a") as file:
                file.write(selected_data_2 + "\n")
            with open(file_names[2],"a") as file:
                file.write(selected_data_3 + "\n")



