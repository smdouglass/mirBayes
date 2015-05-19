%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%


%  This function creates a MATLAB variable that has the genome sequence
%  defined by the user input FASTA file, genomeFileName.  This variable is
%  needed to to compute both entropy and miRNA* presence by future
%  functions, and is required to translate the headers in the alignment
%  file to chromosome numbers.  The genome variable has two fields, Header,
%  which is the FASTA header in the genome file without the '>' character,
%  and Sequence, which contains the sequence for that chromosome.  This
%  variable is returned by the function.

%First get approximate sizes:
function [genome] = readInGenome(genomeFileName)
fid = fopen(genomeFileName, 'r');
if (fid == -1)
    error('Input argument must be the name of a genome sequence in FASTA format')
end
currentIndex = 0;
lengths = [];

while (~feof(fid))
    currentLine = fgetl(fid);
    if (strcmp(currentLine(1), '>'))
        currentIndex = currentIndex + 1;
        lengths(currentIndex) = 0;
    else
        lengths(currentIndex) = lengths(currentIndex) + 1;
    end
end

fclose(fid);
%Done finding lengths of chromosomes for preallocating

%Now actually read in the sequences:
fid = fopen(genomeFileName, 'r');
genome = [];
currentIndex = 0;

while (~feof(fid))
    currentLine = fgetl(fid);
    if (strcmp(currentLine(1), '>'))
        ['done with ' num2str(currentIndex)]
        currentIndex = currentIndex + 1;
        genome(currentIndex).Header = currentLine(2:length(currentLine));
        genome(currentIndex).Sequence = repmat(char(0), 1, lengths(currentIndex)*80);
        currentPosition = 1;
    else
        genome(currentIndex).Sequence(currentPosition:currentPosition+length(currentLine)-1) = currentLine;
        currentPosition = currentPosition + length(currentLine);
    end
end

fclose(fid);