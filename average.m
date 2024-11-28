
%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright �  2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% average.m     This function calculate the average of Vmax for each generation. 

function ave2 = average(pop,popSize)
global VmaxNum;
ave2=mean(pop,2)
ave2(1)=0;
ave2=reshape(ave2,[1,size(ave2,1)]);  