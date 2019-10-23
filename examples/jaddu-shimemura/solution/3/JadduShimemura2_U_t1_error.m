function E = JadduShimemura2_U_t1_error(t1,y2_t1)

U1 = JadduShimemura2_U_t1(t1,t1,y2_t1);
U2 = JadduShimemura2_U_t2(t1,t1,y2_t1);

E = U1-U2;

end