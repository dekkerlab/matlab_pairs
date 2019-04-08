function[upsideDownMatrix] = makeUpsideDownMatrix(originalMatrix)
   [nrows, ncols] = size(originalMatrix);
   upsideDownMatrix = zeros(nrows, ncols);
   for i = 1:nrows
      upsideDownMatrix(i,1:ncols) = originalMatrix(nrows - i + 1, 1:ncols);
   end
end