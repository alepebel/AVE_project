function S = getstim(compfreqs,amps,eventdur,rampdur,sf)

S = zeros(1,eventdur(1)*sf);

% build the sounds
for jj = 1:length(compfreqs)
    S = S + cosramp(amps(jj)*sin(2*pi*compfreqs(jj)*(1/sf:1/sf:eventdur(1))),rampdur,sf);
end
S = S * (1/length(compfreqs));
