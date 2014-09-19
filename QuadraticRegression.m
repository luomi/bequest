function reg = QuadraticRegression(Data,R)


% Align Data
if size(Data,2)>size(Data,1)
    Data=Data';
end

% DegreeWorks

% Function Parameters
NData = size(Data,1);
NVars = size(Data,2);
Coeff_Size=1+NVars+(nchoosek(NVars,2)+NVars);


Regressor_Store=cell(Coeff_Size,1);



%Generate Regressors:
Regressor_Store{1}=ones(NData,1);
count=1;
for i=1:NVars;
    Regressor_Store{i+1}=Data(:,i);
    
    for j=i:NVars
        Regressor_Store{(NVars)+1+count}=   Data(:,i).*Data(:,j);
        count=count+1;
    end
end

X=zeros(NData, Coeff_Size);
for k=1:Coeff_Size
    X(:,k)=Regressor_Store{k};
end



COEFF=(X'*X)\(X'*R);






Constant=COEFF(1);
Linear=COEFF(2:(NVars+1));

Squared=zeros(NVars,NVars);

count=NVars+2;
for i=1:NVars
    for j=i:NVars
        Squared(i,j)=COEFF(count);
        if i~=j
            Squared(i,j)=.5*Squared(i,j);
            Squared(j,i)=Squared(i,j);
        end
        count=count+1;
    end
end

Squared=2.*Squared;



reg = struct('Constant', Constant, 'Linear', Linear, 'Squared', Squared);
end

