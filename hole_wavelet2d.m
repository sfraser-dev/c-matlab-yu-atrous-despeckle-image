%
%---------------hole_wavelet2d.cmex---------------
%
% This uses a single pass 2d kernel based on the
% kronecker product to prevent aliasing. This is
% an improvement upon FAT2_LM which assumes (wrongly)
% that performing 1d scaling upon the rows then the 
% columns does not introduce too much aliasing.
%
% 2-D a trous wavelet scaling (length 3) with zero padding.
%
% Usage: matrix_out=hole_wavelet2d(matrix_in,scl);
%
% Where	'matrix_in' is the input matrix.
%	'scl' is the scale level.
%	'matrix_out' is the processed matrix.
%
% This hole_wavelet2d function can be used to implement the a 
% trous wavelet transformas shown below. 
%
% e.g.
%	% A Trous scaling
%	at1 = hole_wavelet2d(input_image,0);
%	at2 = hole_wavelet2d(at1,1);
%	at3 = hole_wavelet2d(at2,2);
%	at4 = hole_wavelet2d(at3,3);	% Residual image
%
%	% Wavelet levels
%	w1 = input_image - at1;
%	w2 = at1 - at2; 
% 	w3 = at2 - at3;
% 	w4 = at3 - at4;
%
% 	% Perfect reconstruction of input_image
%	reconstruct =  at4 + w4 + w3 + w2 + w1;
%
