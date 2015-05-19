%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%

%  This function takes in the name of the small RNA sequence reads, known
%  miRNAs, and small RNA sequence coverage variables generated by earlier
%  functions.  The user must provide both the variable name and saved file
%  name for each variable.  The function computes the entropy for each
%  mapped position of each small RNA sequence and known miRNA, and returns
%  modified variables for both reads and known miRNAs that include this
%  information as a vector with length equal to the number of mapped
%  positions of each sequence.

function [reads, known_miRNAs] = computeShannonEntropy(sequenceFileName, sequenceVariableName, miRNAFileName, miRNAVariableName, coverageFileName, coverageVariableName)
reads = load(sequenceFileName, sequenceVariableName);
reads = reads.(sequenceVariableName);
known_miRNAs = load(miRNAFileName, miRNAVariableName);
known_miRNAs = known_miRNAs.(miRNAVariableName);
coverage = load(coverageFileName, coverageVariableName);
coverage = coverage.(coverageVariableName);

for i=1:length(reads)
    reads(i).entropy = zeros(1, length(reads(i).positions));
    for j=1:length(reads(i).positions)
        windowStart = reads(i).positions(j) - 50;
        windowEnd = reads(i).positions(j) + 50;
        if (windowStart < 1)
            windowStart = 1;
        end
        if (windowEnd > length(coverage(reads(i).Xsomes(j)).coverage))
            windowEnd = length(coverage(reads(i).Xsomes(j)).coverage);
        end
        windowCoverage = zeros(1, windowEnd-windowStart+1);
        for k=windowStart:windowEnd
            windowCoverage(k-windowStart+1) = coverage(reads(i).Xsomes(j)).coverage(k);
        end
        entropy = 0;
        for k=windowStart:windowEnd
            fractional_occupancy = coverage(reads(i).Xsomes(j)).coverage(k)/sum(windowCoverage);
            if (fractional_occupancy == 0 || log(fractional_occupancy) == 0)
                entropy = entropy + 0;%MATLAB computes NaN, we want 0
            else
                entropy = entropy + (fractional_occupancy*log(fractional_occupancy));
            end
        end
        if (sum(windowCoverage) == 0)
            reads(i).entropy(j) = NaN;%Some known miRNAs may have zero counts, so they have undefined entropy
        else
            reads(i).entropy(j) = (entropy/length(windowCoverage)) * -1;
        end
    end
end

for i=1:length(known_miRNAs)
    known_miRNAs(i).entropy = zeros(1, length(known_miRNAs(i).positions));
    for j=1:length(known_miRNAs(i).positions)
        windowStart = known_miRNAs(i).positions(j) - 50;
        windowEnd = known_miRNAs(i).positions(j) + 50;
        if (windowStart < 1)
            windowStart = 1;
        end
        if (windowEnd > length(coverage(known_miRNAs(i).Xsomes(j)).coverage))
            windowEnd = length(coverage(known_miRNAs(i).Xsomes(j)).coverage);
        end
        windowCoverage = zeros(1, windowEnd-windowStart+1);
        for k=windowStart:windowEnd
            windowCoverage(k-windowStart+1) = coverage(known_miRNAs(i).Xsomes(j)).coverage(k);
        end
        entropy = 0;
        for k=windowStart:windowEnd
            fractional_occupancy = coverage(known_miRNAs(i).Xsomes(j)).coverage(k)/sum(windowCoverage);
            if (fractional_occupancy == 0 || log(fractional_occupancy) == 0)
                entropy = entropy + 0;%MATLAB computes NaN, we want 0
            else
                entropy = entropy + (fractional_occupancy*log(fractional_occupancy));
            end
        end
        if (sum(windowCoverage) == 0)
            known_miRNAs(i).entropy(j) = NaN;%Some known miRNAs may have zero counts, so they have undefined entropy
        else
            known_miRNAs(i).entropy(j) = (entropy/length(windowCoverage)) * -1;
        end
    end
end