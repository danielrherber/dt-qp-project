function E = JadduShimemura1_U_t2_error(t1,t2,y1_t1)

U2 = JadduShimemura1_U_t2(t2);
U3 = JadduShimemura1_U_t3(t2,t1,t2,y1_t1);

E = U2-U3;

end