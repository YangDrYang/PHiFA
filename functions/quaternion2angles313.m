function varargout = quaternion2angles313(varargin)

varargout = cell(1,length(varargin));
for i = 1:length(varargin)

    angles = zeros(3, size(varargin{i},2));
    arates = zeros(3, size(varargin{i},2));
    for j = 1:size(varargin{i},2)
        angles(:,j) = (Q2Eul(norm_quat(varargin{i}(8:11,j)))).*180/pi;%radian -> degree
        arates(:,j) = varargin{i}(12:14,j).*180/pi;%radian/sec -> degree/sec   
    end
    varargout{i} = zeros(13,size(varargin{i},2));
    varargout{i}(1:7,:) = varargin{i}(1:7,:);
    varargout{i}(8:end,:) = [angles;arates];
end