

clear
load('C:\Users\James\Documents\GitHub\may26_onBig2\WNstim_response_OnParasol_RGC.mat')
innerRetinaAll = innerRetina;
% clear innerRetina

load('C:\Users\James\Documents\GitHub\may26_offBig2\WNstim_response_OffParasol_RGC.mat')
innerRetinaAll.mosaic{2,1}= innerRetina.mosaic{1};

load('C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_RGC.mat')
innerRetinaAll.mosaic{3,1}= innerRetina.mosaic{1};

load('C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_RGC.mat')
innerRetinaAll.mosaic{4,1}= innerRetina.mosaic{1};

clear innerRetina
innerRetina = innerRetinaAll;