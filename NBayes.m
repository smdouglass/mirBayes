%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%  The function takes in the small RNA sequence file and variable names and
%  the known miRNA file and variable names to perform Naive Bayes
%  Classification to compute the probability that each small RNA sequence
%  read and known miRNA is a true miRNA.  The user also provides a two
%  column matrix that indicates the bins that counts should be broken up
%  into, with the first column indicating the lower bound of the bin and
%  the second column indicating the upper bound.  For example, values of
%  1,3 would indicate that read counts of 1, 2, and 3 would all be
%  considered indentically for the purposes of classification.  This is to
%  prevent count sparseness from biasing the results.  The final bin should
%  be in the format X,0 where X is the lower bound for the last bin.  This
%  will put all values greater than or equal to X in the last bin.  The
%  user must also provide a vector of enropy break points to indicate where
%  entropy bins should be separated.  This function uses pseudocounts to
%  prevent any variable from giving definate probabilities due to low
%  counts in some variable state(s).  The function returns modified small
%  RNA sequence read and known miRNA variables that have additional 'post'
%  fields.  This field is a vector of equal length as the number of
%  locations a sequence maps to and indicates the probability each mapped
%  location has to being a miRNA under the Naive Bayesian assumptions
%  assuming a uniform prior probability.


function [reads, known_miRNAs] = NBayes(sequenceFileName, sequenceVariableName, knownFileName, knownVariableName, countBins, entropyBreaks)
%entropyBreak is form:  breaks = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 1]; we don't expect it to start with zero, only start with zero if you want to include zero as its own category
reads = load(sequenceFileName, sequenceVariableName);
reads = reads.(genomeVariableName);
known_miRNAs = load(knownFileName, knownVariableName);
known_miRNAs = known_miRNAs.(knownVariableName);

%Counts:
counts = ones(2, length(countBins));%pseudocount; knowns, unknowns; each bin
for i=1:length(known_miRNAs)
    for j=1:length(countBins(:, 1))
        if (known_miRNAs(i).count >= countBins(j, 1) && known_miRNAs(i).count <= countBins(j, 2))
            counts(1, j) = counts(1, j) + 1;
        end
    end
end
for i=1:length(reads)
    for j=1:length(countBins(:, 1))
        if (reads(i).count >= countBins(j, 1) && reads(i).count <= countBins(j, 2))
            counts(2, j) = counts(2, j) + 1;
        end
    end
end

totalReads(1) = sum(counts(1, :));
totalReads(2) = sum(counts(2, :));
for i=1:length(countBins)
   counts(1, i)  = counts(1, i)/totalReads(1);
   counts(2, i)  = counts(2, i)/totalReads(2);
end

%**LENGTH**:
MINIMUM_LENGTH = Inf;
MAXIMUM_LENGTH = 0;
for i=1:length(known_miRNAs)
    if (known_miRNAs(i).length > MAXIMUM_LENGTH)
        MAXIMUM_LENGTH = known_miRNAs(i).length;
    end
    if (known_miRNAs(i).length < MINIMUM_LENGTH)
        MINIMUM_LENGTH = known_miRNAs(i).length;
    end
end
for i=1:length(reads)
    if (reads(i).length > MAXIMUM_LENGTH)
        MAXIMUM_LENGTH = reads(i).length;
    end
    if (reads(i).length < MINIMUM_LENGTH)
        MINIMUM_LENGTH = reads(i).length;
    end
end

lengths = ones(2, MAXIMUM_LENGTH-MINIMUM_LENGTH+1);
for i=1:length(known_miRNAs)
    lengths(1, known_miRNAs(i).length-MINIMUM_LENGTH+1) = lengths(1, known_miRNAs(i).length-MINIMUM_LENGTH+1) + 1;
end
for i=1:length(reads)
    lengths(2, reads(i).length-MINIMUM_LENGTH+1) = lengths(2, reads(i).length-MINIMUM_LENGTH+1) + 1;
end
lengths(1, :) = lengths(1, :)/sum(lengths(1, :));
lengths(2, :) = lengths(2, :)/sum(lengths(2, :));

%**Multiplicity**:
MAXIMUM_MULTIPLICITY = 0;
for i=1:length(known_miRNAs)
    if (length(known_miRNAs(i).positions) > MAXIMUM_MULTIPLICITY)
        MAXIMUM_MULTIPLICITY = length(known_miRNAs(i).positions);
    end
end
for i=1:length(reads)
    if (length(reads(i).positions) > MAXIMUM_MULTIPLICITY)
        MAXIMUM_MULTIPLICITY = length(reads(i).positions);
    end
end

multiplicity = ones(2, MAXIMUM_MULTIPLICITY);
for i=1:length(known_miRNAs)
    multiplicity(1, length(known_miRNAs(i).positions)) = multiplicity(1, length(known_miRNAs(i).positions)) + 1;
end
for i=1:length(reads)
    multiplicity(2, length(reads(i).positions)) = multiplicity(2, length(reads(i).positions)) + 1;
end
multiplicity(1, :) = multiplicity(1, :)/sum(multiplicity(1, :));
multiplicity(2, :) = multiplicity(2, :)/sum(multiplicity(2, :));

%**Entropy**:
counter = 0;
read_entropies = zeros(1, length(reads)*MAXIMUM_MULTIPLICITY);
for i=1:length(reads)
    for j=1:length(reads(i).entropy)
        if (~(reads(i).count == 1 && reads(i).entropy(j) == 0))
            counter = counter + 1;
            read_entropies(counter) = reads(i).entropy(j);
        end
    end
end
read_entropies = read_entropies(1:counter);

counter = 0;
for i=1:length(known_miRNAs)
    for j=1:length(known_miRNAs(i).entropy)
        if (known_miRNAs(i).count ~= 0 && ~(known_miRNAs(i).count == 1 && known_miRNAs(i).entropy(j) == 0))
            counter = counter + 1;
            known_entropies(counter) = known_miRNAs(i).entropy(j);
        end
    end
end

entropies = ones(2, length(entropyBreaks));%knowns, unknowns
entropyBreaks = [-Inf entropyBreaks];
for i=2:length(entropyBreaks)
    for j=1:length(known_entropies)
        if (known_entropies(j) > entropyBreaks(i-1) && known_entropies(j) <= entropyBreaks(i))
            entropies(1, i-1) = entropies(1, i-1) + 1;
        end
    end
    for j=1:length(read_entropies)
        if (read_entropies(j) > entropyBreaks(i-1) && read_entropies(j) <= entropyBreaks(i))
            entropies(2, i-1) = entropies(2, i-1) + 1;
        end
    end
end
entropies(1, :) = entropies(1, :)/sum(entropies(1, :));
entropies(2, :) = entropies(2, :)/sum(entropies(2, :));

%**STAR**:
stars = [0 0];
totalStars = [1 1];%add in pseudocounts to begin with
totalTests = [1 1];
for i=1:length(known_miRNAs)
    if (sum(known_miRNAs(i).star) > 0)
        totalStars(1) = totalStars(1) + 1;
    end
    totalTests(1) = totalTests(1) + 1;
end
for i=1:length(reads)
    if (sum(reads(i).star) > 0)
        totalStars(2) = totalStars(2) + 1;
    end
    totalTests(2) = totalTests(2) + 1;
end
stars(1) = totalStars(1)/totalTests(1);
stars(2) = totalStars(2)/totalTests(2);

%Naive Bayes is: prob(C)*prob(F1|C)*prob(F2|C)*prob(Fn|C):
for i=1:length(reads)
    for j=1:length(reads(i).entropy)
        reads(i).logp(j) = 0;
        reads(i).post(j) = 0;
    end
end

for i=1:length(reads)
    for j=1:length(reads(i).entropy)
        currentProbMir = lengths(1, reads(i).length-MINIMUM_LENGTH+1)*multiplicity(1, length(reads(i).positions));%Do lengths and multiplicity here
        currentProbNotMiR = lengths(2, reads(i).length-MINIMUM_LENGTH+1)*multiplicity(2, length(reads(i).positions));
        if (reads(i).star(j) == 1)
            currentProbMir = currentProbMir * stars(1);%Do miRNA* here
            currentProbNotMiR = currentProbNotMiR * stars(2);
        else
            currentProbMir = currentProbMir * (1 - stars(1));
            currentProbNotMiR = currentProbNotMiR * (1 - stars(2));
        end
        %Do entropy:
        for k=2:length(entropyBreaks)
            if (reads(i).entropy(j) > entropyBreaks(k-1) && reads(i).entropy(j) <= entropyBreaks(k))
                currentProbMir = currentProbMir * entropies(1, k);
                currentProbNotMiR = currentProbNotMiR * entropies(2, k);
            end
        end
        %Do counts:
        for k=1:length(countBins)
            if (reads(i).count >= countBins(k, 1) && reads(i).count <= countBins(k, 2))
                currentProbMir = currentProbMir * counts(1, k);
                currentProbNotMiR = currentProbNotMiR * counts(2, k);
            end
        end
        reads(i).post(j) = currentProbMir/(currentProbMir+currentProbNotMiR);
    end
end

%Knowns:
for i=1:length(known_miRNAs)
    for j=1:length(known_miRNAs(i).entropy)
        currentProbMir = lengths(1, known_miRNAs(i).length-MINIMUM_LENGTH+1)*multiplicity(1, length(known_miRNAs(i).positions));%Do lengths and multiplicity here
        currentProbNotMiR = lengths(2, known_miRNAs(i).length-MINIMUM_LENGTH+1)*multiplicity(2, length(known_miRNAs(i).positions));
        if (known_miRNAs(i).star(j) == 1)
            currentProbMir = currentProbMir * stars(1);%Do miRNA* here
            currentProbNotMiR = currentProbNotMiR * stars(2);
        else
            currentProbMir = currentProbMir * (1 - stars(1));
            currentProbNotMiR = currentProbNotMiR * (1 - stars(2));
        end
        %Do entropy:
        for k=2:length(entropyBreaks)
            if (known_miRNAs(i).entropy(j) > entropyBreaks(k-1) && known_miRNAs(i).entropy(j) <= entropyBreaks(k))
                currentProbMir = currentProbMir * entropies(1, k);
                currentProbNotMiR = currentProbNotMiR * entropies(2, k);
            end
        end
        %Do counts:
        for k=1:length(countBins)
            if (known_miRNAs(i).count >= countBins(k, 1) && known_miRNAs(i).count <= countBins(k, 2))
                currentProbMir = currentProbMir * counts(1, k);
                currentProbNotMiR = currentProbNotMiR * counts(2, k);
            end
        end
        known_miRNAs(i).post(j) = currentProbMir/(currentProbMir+currentProbNotMiR);
    end
end