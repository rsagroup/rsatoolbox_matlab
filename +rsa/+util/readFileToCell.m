% Reads a file line-by-line into a cell array
% CW 1-2010

function cellOut = readFileToCell(fileAddress)

    fileID = fopen(fileAddress, 'r');

    cellOut = {};

    while 1

        lineFeed = fgetl(fileID);

        if ~ischar(lineFeed), break; end%if

        cellOut = [cellOut; {lineFeed}];

    end%while(1)

    fclose(fileID);

end%function
