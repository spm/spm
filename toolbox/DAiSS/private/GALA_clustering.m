function res = GALA_clustering(lJcov,J1, S, distance, A)   

    Nl = S.Nl;
    Nd = length(lJcov)/Nl;
        
    CJcov = GALA_find_localmin(lJcov,Nl,Nd,A,S.threshold);


    for p=1:Nl
        tCJcov=[];
        tCJcov.clvi=CJcov.clvi(p,:);
        tCJcov.maxs=CJcov.maxs{p};
        mxymcor = 0;
        while mxymcor<S.mincorr
            oldsize = length(tCJcov.maxs);
            cursize = length(tCJcov.maxs);
            for i=1:oldsize
                sxycor = corrcoef(J1(tCJcov.clvi{i}+(p-1)*Nd,:)');
                sxycor(sxycor==1)=NaN;
                sxymcor = nanmean(nanmean(sxycor));
                if length(tCJcov.clvi{i})>S.maxsize || sxymcor<S.mincorr
                    cltd = tCJcov.clvi{i};
    %                 ntd = ceil(length(cltd)/20);
                    ntd = 2;
                    Jcovpi = J1(cltd+(p-1)*Nd,:);
                    d = distance(cltd,cltd);
                    Jcd = ones(length(cltd),length(cltd))-corrcoef(Jcovpi');
                    dd = d+S.distratio1*Jcd;
                    Y = squareform(dd);
                    Z = linkage(Y,'complete');
                    CL = cluster(Z,'maxclust',ntd);
                    clm = CL(cltd==tCJcov.maxs(i));
                    tCJcov.clvi{i} = cltd(CL==clm);
                    cld = setdiff(1:max(CL),clm);
                    for j=cld
                        ind = find(CL==j);
                        tCJcov.clvi{cursize+1} = cltd(ind);
                        dd = distance(cltd(ind),cltd(ind));
                        [nu md] = min(sum(dd));
                        tCJcov.maxs(cursize+1) = cltd(ind(md));
                        cursize = cursize+1;
                    end
                end
            end
            i=1;
            while i<=cursize
                if length(tCJcov.clvi{i})<5
                    tCJcov.clvi(i) = [];
                    tCJcov.maxs(i) = [];
                    cursize = length(tCJcov.maxs);
                else 
                    i=i+1;
                end
            end
            xymcor=[];
            for i=1:cursize
                xycor = corrcoef(J1(tCJcov.clvi{i}+(p-1)*Nd,:)');
                xycor(xycor==1)=NaN;
                xymcor(i) = nanmean(nanmean(xycor));
            end
            mxymcor = min(xymcor);
        end
        for i=1:max(length(CJcov.maxs{p}),cursize)
            if i<=cursize
                CJcov.clvi{p,i}=tCJcov.clvi{i};
                CJcov.maxs{p}(i)=tCJcov.maxs(i);
            else
                CJcov.clvi{p,i}=[];
                CJcov.maxs{p}(i)=[];
            end
        end
    end
    
    maxsa = spm_cat(CJcov.maxs);

    ldist = distance(maxsa,maxsa);
    lmaxs = [];
    smaxs = [];
    clind = [];
    for p=1:Nl
        lmaxs = horzcat(lmaxs,CJcov.maxs{p}+Nd*(p-1));
        smaxs(end+1:end+length(CJcov.maxs{p}),:) =...
            horzcat(repmat(p,[length(CJcov.maxs{p}) 1]),CJcov.maxs{p}');
        clind(end+1:end+length(CJcov.maxs{p}),:) =...
            horzcat(repmat(p,[length(CJcov.maxs{p}) 1]),(1:length(CJcov.maxs{p}))');

    end
    lmaxs = lmaxs';

    Jcor = ones(length(lmaxs),length(lmaxs))-corrcoef(J1(lmaxs,:)');
    dd = ldist+S.distratio2*Jcor;
    
    Y = squareform(dd);    
    Z = linkage(Y,S.linkmeth);
    if isfield(S.cluster,'maxclust')
        CL = cluster(Z,'maxclust',S.cluster.maxclust.maxclustsize);
    elseif isfield(S.cluster,'cutoff')
        CL = cluster(Z,'cutoff',S.cluster.cutoff.cutoffthresh,'criterion','distance');
    end

    % remove multiple maxs for one subject in one cluster
    % and merge appropriate basins
    rs=[];
    for cl=1:max(CL)
       list = smaxs(CL==cl,:);
       ind = find(CL==cl);
       lm = lJcov(lmaxs(ind(1)));
       ind_lm = 1;
       for i=2:size(list,1)
           if list(i,1)==list(i-1,1)
               if lJcov(lmaxs(ind(i)))>lm
                   rs = [rs;ind(ind_lm)];
                   CJcov.clvi{clind(ind(i),1),clind(ind(i),2)}=...
                       horzcat(CJcov.clvi{clind(ind(i),1),clind(ind(i),2)},...
                       CJcov.clvi{clind(ind(ind_lm),1),clind(ind(ind_lm),2)});
                   CJcov.clvi{clind(ind(ind_lm),1),clind(ind(ind_lm),2)}=[];
                   lm = lJcov(lmaxs(ind(i)));
                   ind_lm = i;
               else
                   rs = [rs;ind(i)];
                   CJcov.clvi{clind(ind(ind_lm),1),clind(ind(ind_lm),2)}=...
                       horzcat(CJcov.clvi{clind(ind(i),1),clind(ind(i),2)},...
                       CJcov.clvi{clind(ind(ind_lm),1),clind(ind(ind_lm),2)});
                   CJcov.clvi{clind(ind(i),1),clind(ind(i),2)}=[];
               end
           else 
               lm = lJcov(lmaxs(ind(i)));
               ind_lm = i;
           end
       end
    end
    
    
    CL(rs)=[];
    lmaxs(rs)=[];
    smaxs(rs,:)=[];
    clind(rs,:)=[];
    ldist(rs,:) = [];
    ldist(:,rs) = [];
    
    
    CLi = unique(CL);
    
thr = S.Nl-S.similarity;
    
dd=[];list=[];ncl=0;pclvi=[];pmaxs=[];tc=0;
for cl=1:length(CLi)
    if length(lmaxs(CL==CLi(cl)))>=thr
        ncl=ncl+1;
        ind = find(CL==CLi(cl));
        dd = ldist(ind,ind);
        [nu, md] = min(sum(dd));
        list{CLi(cl)} = smaxs(ind,:);
        mvd = list{CLi(cl)}(md,2);
        no = find(~ismember(1:Nl,list{CLi(cl)}(:,1)));
        iclvi = CJcov.clvi{clind(ind(1),1),clind(ind(1),2)};
        uclvi = CJcov.clvi{clind(ind(1),1),clind(ind(1),2)};
        wi = 1;
        while isempty(iclvi)
            wi=wi+1;
            iclvi = CJcov.clvi{clind(ind(wi),1),clind(ind(wi),2)};
            no = sort([no clind(ind(wi-1))]);
            ind(wi-1)=[];
        end
        pclvi{clind(ind(1),1),ncl} = iclvi;
        pmaxs{clind(ind(1),1)}(ncl) = smaxs(ind(1),2);
        for i=2:length(ind)
            pclvi{clind(ind(i),1),ncl} = CJcov.clvi{clind(ind(i),1),clind(ind(i),2)};
            pmaxs{clind(ind(i),1)}(ncl) = smaxs(ind(i),2);
            if isempty(pclvi{clind(ind(i),1),ncl})
                pclvi{clind(ind(i),1),ncl} = iclvi;
                pmaxs{clind(ind(i),1)}(ncl) = smaxs(ind(1),2);
            end
            iclvin = intersect(iclvi,pclvi{clind(ind(i),1),ncl});
            uclvin = union(uclvi,pclvi{clind(ind(i),1),ncl});
        end
        oth = setdiff(1:Nl,list{CLi(cl)}(:,1));
        for i=1:length(oth)
            if length(iclvi)>15
                pclvi{oth(i),ncl} = iclvi;
                pmaxs{oth(i)}(ncl) = smaxs(ind(1),2);
            else
                dd = distance(uclvin,uclvin);
                [nu, md] = min(sum(dd));
                pmaxs{oth(i)}(ncl) = uclvin(md);
                pclvi{oth(i),ncl} = find(distance(uclvin(md),:)<=2);
            end
        end
        
    end    
        
end

res.pclvi = pclvi;
res.pmaxs = pmaxs;



%     figure;
%     ahp = [];
%     for p=1:Nl
% 
%         srcs_disp1 = zeros(Nd,1);
%         for cl=1:size(pclvi,2)
%             srcs_disp1(pclvi{p,cl})=1;
%         end
%         
%         ahp(p)=subplot(1,2,p);
% 
%         cla; axis off
% 
%         fig1 = patch('vertices',vert,'faces',face,'FaceVertexCData',srcs_disp1);
% 
%         set(fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
%         shading interp
%         lighting gouraud
%         %camlight
%         zoom off
%         lightangle(0,270);lightangle(270,0),lightangle(90,0),lightangle(0,45),lightangle(0,135);
%         material([.1 .1 .4 .5 .4]);
%         %view(140,15);
%         caxis([0 1])
%         colormap(jet);
%         hold on;
%         scatter3(vert(pmaxs{p},1),vert(pmaxs{p},2),vert(pmaxs{p},3),15,'r','fill');hold on; 
% 
%     end
%     hlink = linkprop(ahp, {'CameraPosition','CameraUpVector'});
%     key = 'graphics_linkprop';
%     setappdata(ahp(1),key,hlink); 



end
