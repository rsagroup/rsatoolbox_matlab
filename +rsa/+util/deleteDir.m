% FJ 10-2014

function deleteDir(userOptions, Models)

    if (userOptions.deleteTMaps_Dir)
        dirToDel= fullfile(userOptions.rootPath, '/Maps');
        rmdir(dirToDel, 's');
    end
    
    if (userOptions.deleteImageData_Dir)
        dirToDel= fullfile(userOptions.rootPath, '/ImageData');
        rmdir(dirToDel, 's');
    end
    
     if (userOptions.deletePerm)
	    dirToDel= fullfile(userOptions.rootPath, '/Maps/',Models.name, '/perm*.stc');
        delete (dirToDel);  
    end
end
