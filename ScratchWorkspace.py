import arcpy, os

def CreateScratchWorkspace(workspace, defname='scratch'):
    defworkspace = arcpy.env.workspace

    # check if the current path maps to geodatabase
    n = len(workspace)
    if n > 4:
        end = workspace[n - 4: n]  # extract last 4 letters
        if end == ".gdb":  # geodatabase
            workspace = os.path.dirname(workspace)

    arcpy.env.workspace = workspace
    workspaces = arcpy.ListWorkspaces("*", "Folder")
    names = [os.path.basename(w) for w in workspaces]
    i = 0
    name = defname
    while name in names:
        name = defname + str(i)
        i += 1
    arcpy.CreateFolder_management(workspace, name)

    scratchworkspace = workspace + '/' + name
    arcpy.env.workspace = defworkspace

    return scratchworkspace