function extract


load('water level (when dry_ bed level)');
runid1=1;
runid2=length(data.Time);
for i=runid1:runid2
    wld_i=num2str(i);
    wld_name=strcat('wld',wld_i);
    wld=data.Val(i,:,:);
    wld=reshape(wld,size(wld,2),size(wld,3));
    save(wld_name,'wld');
    
end


load('water level');
for i=runid1:runid2
    wl_i=num2str(i);
    wl_name=strcat('wl',wl_i);
    wl=data.Val(i,:,:);
    wl=reshape(wl,size(wl,2),size(wl,3));
    save(wl_name,'wl');
    
end


load('bed level in water level points');
for i=runid1:runid2
    h_i=num2str(i);
    h_name=strcat('h',h_i);
    h=data.Val(i,:,:);
    h=reshape(h,size(h,2),size(h,3));
    save(h_name,'h');
    
end


load('depth averaged velocity');
for i=runid1:runid2
    v_i=num2str(i);
    v_name=strcat('v',v_i);
    v=data.Val(i,:,:);
    v=reshape(v,size(v,2),size(v,3));
    save(v_name,'v');
    
end
  

load('water depth');
for i=runid1:runid2
    d_i=num2str(i);
    d_name=strcat('d',d_i);
    d=data.Val(i,:,:);
    d=reshape(d,size(d,2),size(d,3));
    save(d_name,'d');
    
end


load('total transport');
for i=runid1:runid2
    q_i=num2str(i);
    q_name=strcat('q',q_i);
    q=data.Val(i,:,:);
    q=reshape(q,size(q,2),size(q,3));
    save(q_name,'q');
    
end


end

