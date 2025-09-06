%% define the inner product
function y = innerProduct(a, b, dS_field1, dS_field2)
%     load("config.mat")
%     ds = (length_y/resolutionY) * (length_x/resolutionX); % in m
%     y = ds * ((alpha'.*a')*b);
    length_field1 = length(dS_field1);
%     length_field2 = length(dS_field2);
    y1 = a(1:length_field1,:)' * (b(1:length_field1).*dS_field1);
    y2 = a(length_field1+1:end,:)' * (b(length_field1+1:end).*dS_field2);

%     y1 = a(1:length_field1,:)' * (b(1:length_field1));
%     y2 = a(length_field1+1:end,:)' * (b(length_field1+1:end));
    y = y1 + y2;
%     y = a'*b;
end
