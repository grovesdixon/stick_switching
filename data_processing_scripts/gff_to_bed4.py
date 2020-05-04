#!/usr/bin/env python
##gff_to_bed12.py
##written 9/18/18
##Groves Dixon

#import modules
import argparse
from sys import exit


##############################
###### DEFINE FUNCTIONS ######
##############################

def read_gff(inputGff, featureString, idString):
    """Readin the GFF line by line and collect data for those with blockString in third column.
    Record these in a dictionary keyed to the parent ID for each, parsed from column 9 using parentString.
    
    eg:
    gffDict = {transcript1 : [exonObject1, exonObject2] ...}
    """
    print("\nReading in GFF {}...".format(inputGff))
    bedList = []
    skippedList = []
    totalBlocks=0
    with open(inputGff, 'r') as infile:
        for line in infile:

            #skip comment lines and blank lines
            if line[0]=="#" or not line.strip("\n"):
                continue

            #capture lines indicated by lineIndicator (eg "gene" in 3rd column)
            lineString = line
            line=line.strip("\n").split("\t")
            
            #grab gene components
            scaff = line[0]
            start = line[3]
            stop = line[4]
            strand = line[6]
            description = line[8]

            #record information for block line
            if line[2]==featureString:
                totalBlocks += 1

                #do some checks on this block line
                if idString not in description:
                    print("WARNING. The given parsing string, {}, was not found for this line:".format(idString))
                    print(lineString)
                    print("Suggest checking your GFF and make sure if fits with given arguments")
                    print("Skipping this line.")
                    skippedList.append(line)
                    continue

                #otherwise parse name and store
                else:
                    parsedId = description.split(idString)[1].split(';')[0]
                    bedLine = [scaff, start, stop, parsedId]
                    bedList.append(bedLine)




    
    #print results summary
    print("...\nDone reading GFF.")
    print("\tFound {} total lines indicated with '{}' in the third column".format(totalBlocks, featureString))
    if (len(skippedList)>0):
        print('\tWARNING. {} lines lacked the expected ID string in the description.'.format(len(skippedList)))
    return(bedList)





def output_bed(bedList, outputName):
    """Outputs the data in bed format"""
    print("\nWriting out results to {}...".format(outputName))
    written=0
    with open(outputName, 'w') as out:
        for bl in bedList:
            written+=1
            outString = '\t'.join(bl)
            if written==1:
                out.write(outString)
            else:
                out.write('\n' + outString)



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Generate a bed file from a gff for a given feature type.

    Output follows ensembl's description of bed file:
    1. chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
    2. chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
    3. chromEnd - End position of the feature in standard chromosomal coordinates
    4. name - Label to be displayed under the feature, if turned on in "Configure this page".

    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-gff', required = True, dest = 'input_gff', help = 'The the input gff')
    parser.add_argument('-o', required = True, dest = 'output_name', help = 'Name for output bed file')
    parser.add_argument('-feature', required = True, dest = 'feature_string', help = 'String for the type of feature to make bed lines from. Should match third column in gff (gene, mRNA, CDS, etc.)')
    parser.add_argument('-IDstring', required = True, default='ID', dest = 'id_string', help = 'String that indicates the id label in the description column you want to use (default = "ID")')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    inputGff = args.input_gff
    outputName = args.output_name
    featureString = args.feature_string
    idString = args.id_string + "="



    #---- RUN FUNCTIONS ----#
    bedList = read_gff(inputGff, featureString, idString)
    output_bed(bedList, outputName)


