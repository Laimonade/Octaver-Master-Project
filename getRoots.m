



function complx = getRoots(div, count)
    roots = [];
    if isempty(roots)
        roots = complex(zeros(count, 1, 'single'));
        cosInc = cos((2*pi) / div);
        sinInc = sin((2*pi) / div);

        roots(1) = complex(1, 0);
        roots = complex(roots);
        lre = 1.0;
        lim = 0.0;
        for i = 1:count
            re = cosInc + lre - sinInc * lim;
            im = sinInc * lre + cosInc * lim;
            lre = re;
            lim = im;
            roots(i) = complex(single(re), single(im));
        end
    end
    complx = roots;
end


% Cmplx[] getRoots(int div, int count) {
%         Cmplx[] roots = nToRoots.get(div);
%         if (roots == null) {
%             roots = Cmplx.newArray(count);
%             double cosInc = Math.cos(Math.PI * 2.0 / div);
%             double sinInc = Math.sin(Math.PI * 2.0 / div);
%             roots[0].set(1.0f, 0.0f);
%             double lre = 1.0;
%             double lim = 0.0;
%             for (int i = 1; i < count; i++) {
%                 /*
%                  * cos(A+B) = cos(A)cos(B) - sin(A)sin(B)
%                  * sin(A+B) = sin(A)cos(B) + cos(A)sin(B)
%                  * with B = first root and A = previous root
%                  */
%                 double re = cosInc * lre - sinInc * lim;
%                 double im = sinInc * lre + cosInc * lim;
%                 lre = re;
%                 lim = im;
%                 roots[i].set((float) re, (float) im);
%             }
%             nToRoots.put(div, roots);