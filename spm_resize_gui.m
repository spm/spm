
% --- Executes when typeAssignmentPanel is resized.
function spm_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to typeAssignmentPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

children=findall(hObject);

for i=1:length(children)
    if children(i)~=hObject
        try
            set(children(i), 'Position', get(children(i), 'Position'));
        end
    end
end
