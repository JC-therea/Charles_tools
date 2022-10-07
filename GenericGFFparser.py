import re

gff = "Whatever.gff"
output = "WhateverOut.gff"

Raw_file_o = open(gff, "r").readlines()
Clean_file_o = open(output, "w+")

for line in Raw_file_o:
    line_without_n = line.replace("\n", "")
    line_splitted = re.split("\t", line_without_n)

    if len(line_splitted) > 1:
        Chr, Source, Feature, Start, End, Score, Strand, Frame, Attributes = line_splitted

        # Here yo do the stuffs
        if Feature == "CDS" or Feature == "exon":
            if Attributes.find(",") >= 0:
                New_attr = re.sub(";.*", "", Attributes.replace("Parent=", "")).split(",")
                for Exon in New_attr:
                    line_splitted[-1] = "Parent=" + Exon
                    new_line = "\t".join(line_splitted) + "\n"
                    Clean_file_o.writelines(new_line)

            else:
                Clean_file_o.writelines(line)
        else:
            Clean_file_o.writelines(line)

Clean_file_o.close()