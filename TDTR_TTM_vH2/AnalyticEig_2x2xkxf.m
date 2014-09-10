function [Eig, Lambda] = AnalyticEig_2x2xkxf(A)
% AnalyticEig gives the analytic solution to the 2x2 eigenvalue problem
% parallelized across "k" AND "f" spaces. That is, it solves kxf 2x2 eigenvalue
% problems in parallel.

% Each problem is of the form A(:,:,i,j).*Eig = Lambda.*Eig, where ...
% Lambda are the eigenvalues and
% Eig are the normalized eigenvectors.
% Eig(:,1,:,:) is the column eigenvector for Lambda(1,1,:,:),
% Eig(:,2,:,:) is the column eigenvector for Lambda(2,1,:,:).
% (A - Lambda.*I).*Eig = 0
% ThermalMatrix .* x = - Kchannels .* Lambda .* x
% Solution from the Internet:
% https:././bearspace.baylor.edu./Vince_Cronin./www./PBO_ed./FindingEigenvectors.pdf

%% Computes eigenvalues.
Lambda(2,1,:,:) = ((A(1,1,:,:) + A(2,2,:,:)) + sqrt((4.*A(1,2,:,:).*A(2,1,:,:)) + (A(1,1,:,:) - A(2,2,:,:)).^2)) ./ 2;
Lambda(1,1,:,:) = ((A(1,1,:,:) + A(2,2,:,:)) - sqrt((4.*A(1,2,:,:).*A(2,1,:,:)) + (A(1,1,:,:) - A(2,2,:,:)).^2)) ./ 2;

%% Computes eigenvectors. The if/else structures handle degenerate cases.
% I don't know how to parallelize the normalization, so I skipped that
% step.

%Norm1 = sqrt(1+((Lambda(1,1,:,:)-A(1,1,:,:))./A(1,2,:,:)).^2);
%Norm2 = sqrt(1+((Lambda(2,1,:,:)-A(1,1,:,:))./A(1,2,:,:)).^2);

Eig = ones(2,2,length(A(1,1,:,1)),length(A(1,1,1,:)));

% Prevents errors from singleton dimension in Lambda(1,1,:,1).
temp = ones(2,1,length(A(1,1,:,1)),length(A(1,1,1,:)));
temp(:,1,:,:) = Lambda;
Lambda = temp;

if A(1,2,:,:) ~= 0
    Eig(1,1,:,:) = 1;
    Eig(2,1,:,:) = ((Lambda(1,1,:,:) - A(1,1,:,:)) ./ A(1,2,:,:));
else
    if A(2,2,:,:) - Lambda(1,1,:,:) ~= 0
        Eig(1,1,:,:) = 1;
        Eig(2,1,:,:) = - A(2,1,:,:) ./ (A(2,2,:,:) - Lambda(1,1,:,:));
    else
        Eig(1,1,:,:) = 0; 
        Eig(2,1,:,:) = 1;
    end
end
%[absnorm,index] = max(abs(Eig(:,1,:,:)));
%Norm = absnorm.*sign(Eig(index,1,:,:));
%Eig(:,1,:,:) = Eig(:,1,:,:)./Norm;

if A(1,2,:,:) ~= 0
    Eig(1,2,:,:) = 1;
    Eig(2,2,:,:) = ((Lambda(2,1,:,:) - A(1,1,:,:)) ./ A(1,2,:,:));
else
    if A(2,2,:,:) - Lambda(2,1,:,:) ~= 0
        Eig(1,2,:,:) = 1;
        Eig(2,2,:,:) = - A(2,1,:,:) ./ (A(2,2,:,:) - Lambda(2,1,:,:));
    else
        if Eig(1,1,:,:) == 0 % Eig1 is a multiple of (0,1,:,:)
            Eig(1,2,:,:) = 0; % the only orthogonal vector to (0,1,:,:) is (1,0)
            Eig(2,2,:,:) = 1;
        else % Eig1 is not a multiple of (0,1,:,:)
            Eig(1,2,:,:) = 0; % (0,1,:,:) is orthogonal to non-multiples of (0,1,:,:)
            Eig(2,2,:,:) = 1;
        end
    end
end
%[absnorm,index] = max(abs(Eig(:,2,:,:)));
%Norm(:) = absnorm.*sign(Eig(index,2,:,:));
%Eig(:,2,:,:) = Eig(:,2,:,:)./Norm;

%if Eig(1,1,:,:) == Eig(1,2,:,:) && Eig(2,1,:,:) == Eig(2,2,:,:)
%    Eig(:,2,:,:) = -Eig(:,1,:,:); % degenerate matrix; may as well make the eigenvectors _look_ different.
%end

end

