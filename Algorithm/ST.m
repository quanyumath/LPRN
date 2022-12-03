function [T1SQ, I1, T2SQ, I2, T3SQ, I3] = ST(X)
[n1, n2, n3] = size(X);
diagT_3 = eye(n3, n3-1);
for i = n3:-1:2
    diagT_3(i, :) = diagT_3(i, :) - diagT_3(i-1, :);
end
T3(1, :, :) = diagT_3;
for i = 2:n1
    T3(i, :, :) = zeros(size(diagT_3));
end
T3_bar = fft(T3, [], 1);
for u = 1:ceil((n1 + 1)/2)
    TT3 = squeeze(T3_bar(u, :, :));
    T3SQ{u} = TT3 * TT3';
end
I3 = eye(size(T3SQ{1}));
%-----------------------------------
diagT_2 = diag(ones(1, n2)) + diag(repmat(-1/2, 1, n2-1), 1) + diag(repmat(-1/2, 1, n2-1), -1);
diagT_2(1, 2) = -1;
diagT_2(n2, n2-1) = -1;
% diagT_2=eye(n2,n2-1);
% for i=n2:-1:2
%     diagT_2(i,:)=diagT_2(i,:)-diagT_2(i-1,:);
% end
T2(:, :, 1) = diagT_2';
for i = 2:n3
    T2(:, :, i) = zeros(size(diagT_2));
end
T2_bar = fft(T2, [], 3);
for u = 1:ceil((n3 + 1)/2)
    TT2 = squeeze(T2_bar(:, :, u));
    T2SQ{u} = TT2 * TT2';
end
I2 = eye(size(T2SQ{1}));
%-----------------------------------
diagT_1 = diag(ones(1, n1)) + diag(repmat(-1/2, 1, n1-1), 1) + diag(repmat(-1/2, 1, n1-1), -1);
diagT_1(1, 2) = -1;
diagT_1(n1, n1-1) = -1;
% diagT_1=eye(n1,n1-1);
% for i=n1:-1:2
%     diagT_1(i,:)=diagT_1(i,:)-diagT_1(i-1,:);
% end
T1(:, 1, :) = diagT_1;
for i = 2:n2
    T1(:, i, :) = zeros(size(diagT_1'));
end
T1_bar = fft(T1, [], 2);
for u = 1:ceil((n2 + 1)/2)
    TT1 = squeeze(T1_bar(:, u, :));
    T1SQ{u} = TT1' * TT1;
end
I1 = eye(size(T1SQ{1}));

end