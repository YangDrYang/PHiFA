function diffAng = ephAngdiff(varargin)

nEpoch = length(varargin{1}(1,:));
diffAng = zeros(nEpoch,3);
for i = 1:nEpoch
    diffAng(i,1) = angdiff(varargin{1}(8,i)/180*pi, varargin{2}(8,i)/180*pi)/pi*180;
    diffAng(i,2) = angdiff(varargin{1}(9,i)/180*pi, varargin{2}(9,i)/180*pi)/pi*180;
    diffAng(i,3) = angdiff(varargin{1}(10,i)/180*pi, varargin{2}(10,i)/180*pi)/pi*180;
end
