clear all; clc;
disp('start')
[Iw2007,Rw2007,P2007,V2007,Y2007,W2007,USDRate2007,GDPError2007, StaticModelError2007] = StaticModelConstruction(2007);
n = size(Iw2007,1)-1;
[Ia2007,Ra2007] = aggregation(Iw2007,Rw2007,[ones(n,1); 2]);
[Ir2007,Rr2007] = reduction(Iw2007,Rw2007,[20; 35]);
Vb = [1321903; 2012151; 2358546; 2698867; 3512183; 4284803; 5983036];
Wb = [258500; 422600; 490500; 637900; 817000; 1057900; 1327600];
Vb2007 = Vb(end)*USDRate2007;
Wb2007 = Wb(end)*USDRate2007;
rf = 0.1;
[rg,Cn,Rc,rw,Rh,Prb,Pcb,Rntb,Rsb] = IndexesComputation(Iw2007,V2007, Y2007,W2007,Vb2007,Wb2007,rf);
[Prr,Rsp,Rnt] = EfficiencyAnalysis(Ra2007,Ia2007,W2007);
Iw = zeros(n+1,7);
Rw = zeros(n+1,n+1,7);
Dw = zeros(n+1,n+1,6);
Ia = zeros(2,7);
Ra = zeros(2,2,7);
Da = zeros(2,2,6);
Ir = zeros(7,7);
Rr = zeros(7,7,7);
Dr = zeros(7,7,6);
Dl = zeros(7,7,6);

for i = 1:7
    disp('etap')
    disp(i)
[Iw(:,i),Rw(:,:,i),~,~,~,~,~,~,~] = StaticModelConstruction(2000+i);
[Ia(:,i),Ra(:,:,i)] = aggregation(Iw(:,i),Rw(:,:,i),[ones(n,1); 2]);
[Ir(:,i),Rr(:,:,i)] = reduction(Iw(:,i),Rw(:,:,i),[2; 8; 18; 20;...
30; 31; 35]);
if i ~= 1
Dw(:,:,i-1) = DifferentialModelConstruction(Iw(:,i), Iw(:,i-1),Rw(:,:,i));
Da(:,:,i-1) = DifferentialModelConstruction(Ia(:,i), Ia(:,i-1),Ra(:,:,i));
Dr(:,:,i-1) = DifferentialModelConstruction(Ir(:,i), Ir(:,i-1),Rr(:,:,i));
Dl(:,:,i-1) = DifferenceModelConstruction(Ir(:,i), Ir(:,i-1),Rr(:,:,i));
end
end

[pIw,pIwError,pIwErrorMean] = DifferentialModelPrognosis(Iw(:,2:end),Dw);
disp(pIwError)
disp(pIwErrorMean)
[pIa,pIaError,pIaErrorMean] = DifferentialModelPrognosis(Ia(:,2:end),Da);
disp(pIaError)
disp(pIaErrorMean)
[pIr,pIrError,pIrErrorMean] = DifferentialModelPrognosis(Ir(:,2:end),Dr);
disp(pIrError)
disp(pIrErrorMean)
[pIl,pIlError,pIlErrorMean] = DifferenceModelPrognosis(Ir(:,2:end),Dl);
disp(pIlError)
disp(pIlErrorMean)

function [pI,pError,pErrorMean] = DifferenceModelPrognosis(I,D)
[n,m] = size(I);
pI = zeros(n,m-1);
pError = zeros(n,m-1);
for i = 1:m-1
pI(:,i) = D(:,:,i)*I(:,i);
pError(:,i) = abs(pI(:,i)-I(:,i+1))./I(:,i+1);
end
pErrorMean = mean(pError,2);
end

function [pI,pError,pErrorMean] = DifferentialModelPrognosis(I,D)
[n,m] = size(I);
pI = zeros(n,m-1);
pError = zeros(n,m-1);
for i = 1:m-1
pI(:,i) = expm(D(:,:,i))*I(:,i);
pError(:,i) = abs(pI(:,i)-I(:,i+1))./I(:,i+1);
end
pErrorMean = mean(pError,2);
end

function D = DifferenceModelConstruction(I,ILast,R)
M = eye(size(I,1))+diag((I-ILast)./I);
D = M*R;
end

function D = DifferentialModelConstruction(I,ILast,R)
M = diag((I-ILast)./I);
D = M*R;
end

function [Prr,Rsp,Rnt] = EfficiencyAnalysis(Ra,Ia,W)
Rp = Ra(1,1);
Vr = Ra(2,1);
Is = Ia(1,1);
Wr = sum(W)/Is;
Prr = Vr - Wr;
Rsp = Rp + Wr;
Rnt = Prr/Rsp;
end

function [rg,Cn,Rc,rw,Rh,Prb,Pcb,Rntb,Rsb] = IndexesComputation(Iw,V,Y,W,Vb,Wb,rf)
sV = sum(V);
sY = sum(Y);
sW = sum(W);
GDP = Iw(end);
rg = Vb/GDP;
Cn = Vb + sW*(1-rf);
Rc = Cn/sV;
rw = sW/sV;
Rh = 1 - Rc;
Prb = Vb - Wb;
Pcb = sY+Wb;
Rntb = Prb/Pcb;
Rsb = Pcb/GDP;
end

function [Ir,Rr] = reduction(Iw,Rw,unreducedIndustriesNumbers)
n = size(Iw,1)-1;
reducedIndustriesNumbers = setdiff(1:n+1,unreducedIndustriesNumbers);
n1 = length(reducedIndustriesNumbers);
R11 = Rw(unreducedIndustriesNumbers,unreducedIndustriesNumbers);
R12 = Rw(unreducedIndustriesNumbers,reducedIndustriesNumbers);
R21 = Rw(reducedIndustriesNumbers,unreducedIndustriesNumbers);
R22 = Rw(reducedIndustriesNumbers,reducedIndustriesNumbers);
Ir = Iw(unreducedIndustriesNumbers,1);
Rr = R11 + R12*inv(eye(n1)-R22)*R21;
end

function [Ia,Ra] = aggregation(Iw,Rw,aggregateIndustriesNumber)
m = max(aggregateIndustriesNumber);
aggregateIndustry = cell(1,m);
Ia = zeros(m,1);
Ra = zeros(m);
industryAmount = 0;
for i = 1:m
aggregateIndustry{i} = find(aggregateIndustriesNumber==i);
Ia(i,1) = sum(Iw(aggregateIndustry{i},1));
industryAmount = industryAmount+length(aggregateIndustry{i});
end
if industryAmount ~= size(aggregateIndustriesNumber,1)
errordlg('error');
end
for i = 1:m
for j = 1:m
R1 = Rw(aggregateIndustry{i},aggregateIndustry{j});
I1 = Iw(aggregateIndustry{j});
Ra(i,j) = sum(R1*I1)/Ia(j,1);
end
end
end

function [Iw,Rw,P,V,Y,W,USDRate,GDPError,StaticModelError] = StaticModelConstruction(year)
IOTableWorksheet = num2str(year);
CompensationTableRow = char(year-1926);
ExchangeTableRow = char(year-1928);
IO=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'E7:AL40');
ImportsIO=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'E42:AL76');
ValueAdded=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'E78:AL83');
Taxes=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'AN78:AR83');
Consumption=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'AN7:AR40');
Exports=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'AS7:AS40');
Output=xlsread('JPN_NIOT.xlsx', IOTableWorksheet, 'AT7:AT40');
Compensation = xlsread('Compensation.xlsx', 'DATA', [CompensationTableRow '19911:' CompensationTableRow '19944']);
USDRate = xlsread('Exchange.xlsx','EXR',[ExchangeTableRow '27']);
RealJGDP = xlsread('JPN_GDP.xlsx','Data','E2:S2');
P = IO;
V = sum(ValueAdded,1);
Y = sum(Consumption,2);
Tx = sum(sum(Taxes,1),2);
Im = sum(ImportsIO,1);
Ex = Exports;
I = Output;
realGDP = RealJGDP*1000*USDRate;
W = Compensation*USDRate;
GDP = sum(V) + Tx;
GDPError = abs(realGDP-GDP)/realGDP;
n = size(I,1);
Iw = I - Im';
Yw = Y + Ex - Im';
R = zeros(n);
for i =1:n
R(:,i) = P(:,i)/Iw(i,1);
end
Iw = [Iw; GDP];
Vwr = V./Iw(1:n)';
Wwr = V./Iw(1:n)';
Ywr = Yw/Iw(n+1);
Txr = Tx/Iw(n+1);
Rw = [R Ywr; Vwr Txr];
StaticModelError = norm(Iw-Rw*Iw);
end