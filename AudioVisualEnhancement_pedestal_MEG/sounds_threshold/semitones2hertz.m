function f2 = semitones2hertz(f1, semitones)

% f2 = semitones2hertz(f1, semitones)
%
% Inputs:
%	f1        - reference frequency [in Hz], e.g. refFreq = 700;
%	semitones - change from the reference [in semitones]. e.g. semitones = 3;
%
% Output:
%	f2 - new frequency
%
% Description: The programs adds or subtracts the given number of semitones to/from
% a reference frequency (f1).
% -----------------------------------------------------------------------------
% M. Henry, B. Herrmann, Email: bherrmann@cbs.mpg.de, 2011-11-22

nInputs = nargin;
if nInputs < 2, fprintf('Error: Two inputs need to be defined!\n'), f2 = []; return; end;

f2 = f1 * exp(log(2)*semitones/12);

return;