function code=generateCode(type,PRN,Nu,numBand)
if type==1
    load('weil10230_signed.mat');
    code=weil10230_signed(PRN,:);
elseif type==0
    code=generateCAcode(PRN);
else
    load('weil10230_signed.mat');
    code=weil10230_signed(PRN,1:Nu*numBand);
end