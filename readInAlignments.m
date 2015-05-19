%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%


%  This function reads in the Bowtie alignments for both the small RNA
%  sequence reads and known miRNAs that are needed for the Naive Bayes
%  Classification.  The Bowtie alignments need to be in Bowtie's '.bow'
%  format.  The function also requires the file name and variable name for
%  the MATLAB genome variable created by readInGenome.m.  This script reads
%  in the alignments to determine the sequence and counts for each small
%  RNA sequence read, and determines the length and multiplicity from these
%  alignments.  All trusted known miRNA sequences should be aligned for
%  this function, both expressed and non-expressed.  The function will
%  determine the known miRNA counts from the small RNA sequence read
%  alignments.  The function returns variables for both the small RNA
%  sequence reads and known miRNAs that will be used in all forthcoming
%  functions of the Naive Bayes Classifier.

function [reads, known_miRNAs] = readInAlignments(bowtieReadFileName, bowtieMiRNAFileName, genomeFileName, genomeVariableName)
genome = load(genomeFileName, genomeVariableName);
genome = genome.(genomeVariableName);
for i=1:length(genome)
    Xsomemap.(genome(i).Header) = i;
end

%First find the number of unique reads:
numSeqs = 0;
fid = fopen(bowtieReadFileName, 'r');
if (fid == -1)
    error('First input argument must be the name of the Bowtie output file for the small RNA sequence reads')
end

lastID = 'notanid';
while (~feof(fid))
    currentLine = fgetl(fid);
    currentCols = textscan(currentLine, '%s', 'delimiter', '\t', 'MultipleDelimsAsOne', 1);
    currentID = char(currentCols{1}(1));
    if (~strcmp(currentID, lastID))
        numSeqs = numSeqs + 1;
    end
    lastID = currentID;
end
fclose(fid);


fid = fopen(bowtieReadFileName, 'r');
reads(numSeqs).name = '';
reads(numSeqs).positions = [];
reads(numSeqs).Xsomes = [];
reads(numSeqs).length = 0;
reads(numSeqs).count = 0;
readMap = containers.Map;

lastID = 'notanid';
idCounter = 0;
while (~feof(fid))
    currentLine = fgetl(fid);
    currentCols = textscan(currentLine, '%s', 'delimiter', '\t', 'MultipleDelimsAsOne', 1);
    currentSeq = char(currentCols{1}(5));
    if (readMap.isKey(currentSeq))
        if (~strcmp(currentID, lastID))
            reads(readMap(currentSeq)).count = reads(readMap(currentSeq)).count + 1;
        end
    else
        currentXsome = char(currentCols{1}(3));
        currentXsomeIndex = Xsomemap.(currentXsome);
        currentStart = str2double(char(currentCols{1}(4)));
        if (strcmp(currentID, lastID))
            withinCounter = withinCounter + 1;
            reads(idCounter).positions(withinCounter) = currentStart + 1;%+1 to compensate for bowtie's 0-based alignment
            reads(idCounter).Xsomes(withinCounter) = currentXsomeIndex;
        else
            withinCounter = 1;
            idCounter = idCounter + 1;
            readMap(currentSeq) = idCounter;
            currentLength = length(currentSeq);
            reads(idCounter).name = currentSeq;
            reads(idCounter).positions(withinCounter) = currentStart + 1;%+1 to compensate for bowtie's 0-based alignment
            reads(idCounter).Xsomes(withinCounter) = currentXsomeIndex;
            reads(idCounter).length = currentLength;
            reads(idCounter).count = 1;
        end
    end
    lastID = currentID;
end
fclose(fid);

for i=1:length(reads)
    reads(i).multiplicity = length(reads(i).positions);
    reads(i).star = [];
    reads(i).entropy = [];
end

fid = fopen(bowtieMiRNAFileName, 'r');
if (fid == -1)
    error('Second input argument must be the name of the Bowtie output file for the known miRNAs')
end

known_miRNAs = [];
while (~feof(fid))
    currentLine = fgetl(fid);
    currentCols = textscan(currentLine, '%s', 'delimiter', '\t', 'MultipleDelimsAsOne', 1);
    currentSeq = char(currentCols{1}(5));
    currentXsome = char(currentCols{1}(3));
    currentXsomeIndex = Xsomemap.(currentXsome);
    currentStart = str2double(char(currentCols{1}(4)));
    if (strcmp(currentID, lastID))
        withinCounter = withinCounter + 1;
        known_miRNAs(length(known_miRNAs)).positions(withinCounter) = currentStart + 1;%+1 to compensate for bowtie's 0-based alignment
        known_miRNAs(length(known_miRNAs)).Xsomes(withinCounter) = currentXsomeIndex;
    else
        withinCounter = 1;
        currentLength = length(currentSeq);
        known_miRNAs(length(known_miRNAs)).name = currentSeq;
        known_miRNAs(length(known_miRNAs)).positions(withinCounter) = currentStart + 1;%+1 to compensate for bowtie's 0-based alignment
        known_miRNAs(length(known_miRNAs)).Xsomes(withinCounter) = currentXsomeIndex;
        known_miRNAs(length(known_miRNAs)).length = currentLength;
        known_miRNAs(length(known_miRNAs)).count = reads(readMap(currentSeq)).count;
    end
end