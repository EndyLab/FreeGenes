import json
import os
import glob
import re
import csv


def file_list(path):
    counter = 0
    path_number = []
    for files in glob.glob(path):
        counter += 1
        file_name = files.split("/")[-1]
        print("{}. {}".format(counter,file_name))
        path_number.append(files)
    number = input("Which file: ")
    number = int(number) - 1
    return path_number[number]

def genelist_to_list(file_name):
    with open(file_name,"r") as data_file:
        data = data_file.read()
    essential = list(filter(None, data.split("\n")))
    return essential

essential_list = file_list("./../genome/genome_sequences/*.txt")
essential_list = set(genelist_to_list(essential_list))

bioe80_list = set()
for file in glob.glob("./data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    bioe80_list.add(data["info"]["documentation"]["gene_name"])



for gene in bioe80_list^essential_list:
    print(gene)

#with open(name,"w") as write_file:
        #write_file.write(data)
