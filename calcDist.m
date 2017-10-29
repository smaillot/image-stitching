function d = calcDist(H,match1,match2)
    n = length(match1);
    proj = H*[match1' ; ones(1,n)];
    proj = proj(1:2,:)./repmat(proj(3,:),2,1);
    d = sum((match2'-proj).^2,1);

end