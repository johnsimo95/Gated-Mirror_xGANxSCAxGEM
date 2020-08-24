function [port] = ImportCOMSOL_Z(filename)
%IMPORTCOMSOL_Z Imports the port impedance from COMSOL to MATLAB
%   Header variable : f, z, s11, v, i
data=readtable(filename,'HeaderLines',5);
data=table2array(data);

port.freq = data(:,1);
port.z = data(:,2);
port.s11 =  data(:,3);
port.v =  data(:,4);
port.i = data(:,5);

end

