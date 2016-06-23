% Function to read fasta files and get HKY distances etc using MBEToolbox
function dbranch = getMBEdata(name, dtype)

% Assumptions and modifications
% - uses UPGMA tree construction

% Read the fasta file into MBE format
[aln] = readfasta3(name, 1, 0); % for non-coding DNA/RNA

% For each type get the distances between sequences
switch(dtype)
    case 0
        % Jukes-Cantor distance
        [D, ~]= dn_jc(aln);
        disp('Using Jukes-Cantor');
    case 1
        % HKY distance
        [D, ~]= dn_hky(aln);  
        disp('Using HKY');
    case 2
        % Log-Det distance
        D = dn_logdet(aln);
        disp('Using Log-Det');
    case 3
        % Tamura 92 distance
        D = dn_tamura92(aln);
        disp('Using Tamura 92');
    case 4
        % Kimura 80 distance
        [D, ~] = dn_k2p(aln);
        disp('Using Kimura 80');
    case 5
        % Tajima & Nei 84 distance
        [D, ~] = dn_tajima_nei84(aln);
        disp('Using Tajima-Nei 84');
end

% Construct the UPGMA tree and get the branch distances
[topology] = plotupgma(D);
dbranch = topology(:, 4);