% Quick scaling argument

clear
load isoDO_shelf.mat

coeffs=polyfit(2003:2023,isoDO.AmnTemp(2,:),1);
tempRange=polyval(coeffs,[2003 2023]);

salRange=[nanmean(isoDO.SA(2,:)) nanmean(isoDO.SA(2,:))];

[~,sol0]=aou(salRange(1),tempRange(1),isoDO.fittedy(1));
pO20=isoDO.fittedy(1)/sol0;

[~,sol]=aou(salRange(1),tempRange(2),isoDO.fittedy(end));
pO2=isoDO.fittedy(end)/sol;

tempFactor=exp(0.4 / 8.617e-5 * (1 ./ tempRange(2) - 1 ./ tempRange(1)));
