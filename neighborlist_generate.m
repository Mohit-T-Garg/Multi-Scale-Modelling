function neighborlist_generate(x)
    global N;
    global rcut;
    global eps_ver;
    global MAX_NEIGHBOR;
    global x_last;
    global Num;
    global Neighborlist;
    global dim;
    global H;
    global H_mod;

    delta = x - x_last;
    deltamax1 = -500;
    deltamax2 = -500;

    for i = 1:N
        distance = norm(delta(dim * (i - 1) + 1:dim * i));
        if distance > deltamax1
            deltamax2 = deltamax1;
            deltamax1 = distance;
        elseif distance > deltamax2
            deltamax2 = distance;
        end
    end

    if (N == 1)
        deltamax2 = 0;
    end

    if deltamax1 + deltamax2 > eps_ver * rcut
        %disp('neighborlist being updated');
        x_last = x;
        Num = zeros(N, 1);
        Neighborlist = zeros(N, MAX_NEIGHBOR, dim + 1);
        n = zeros(dim, 1);

        if dim == 3
            for alpha = 1:N
                for count1 = -2:2
                    for count2 = -2:2
                        for count3 = -2:2
                            n = [count1; count2; count3];
                            for beta = 1:N
                                if ~(all(n == [0; 0; 0]) && beta == alpha)
                                    imagepos = zeros(dim, 1);

                                    for k = 1:dim
                                        imagepos(k) = x(dim * (beta - 1) + k);

                                        for l = 1:dim
                                            imagepos(k) = imagepos(k) + H_mod(k, l) * n(l);
                                        end
                                    end

                                    dist = norm(x(dim * (alpha - 1) + 1:dim * alpha) - imagepos);

                                    if dist < rcut * (1 + eps_ver)
                                        Num(alpha) = Num(alpha) + 1;
                                        Neighborlist(alpha, Num(alpha), 1) = beta;

                                        for i = 1:dim
                                            Neighborlist(alpha, Num(alpha), i + 1) = n(i);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        if dim == 2
            for alpha = 1:N
                for count1 = -2:2
                    for count2 = -2:2
                        n = [count1; count2];

                        for beta = 1:N
                            if ~(all(n == [0; 0]) && beta == alpha)
                                imagepos = zeros(dim, 1);

                                for k = 1:dim
                                    imagepos(k) = x(dim * (beta - 1) + k);

                                    for l = 1:dim
                                        imagepos(k) = imagepos(k) + H_mod(k, l) * n(l);
                                    end
                                end

                                dist = norm(x(dim * (alpha - 1) + 1:dim * alpha) - imagepos);

                                if dist < rcut * (1 + eps_ver)
                                    Num(alpha) = Num(alpha) + 1;
                                    Neighborlist(alpha, Num(alpha), 1) = beta;

                                    for i = 1:dim
                                        Neighborlist(alpha, Num(alpha), i + 1) = n(i);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        if dim == 1
            for alpha = 1:N
                for count1 = -2:2
                    n = count1;

                    for beta = 1:N
                        if ~(n == 0 && beta == alpha)
                            imagepos = zeros(dim, 1);

                            for k = 1:dim
                                imagepos(k) = x(dim * (beta - 1) + k);

                                for l = 1:dim
                                    imagepos(k) = imagepos(k) + H_mod(k, l) * n(l);
                                end
                            end

                            dist = norm(x(dim * (alpha - 1) + 1:dim * alpha) - imagepos);

                            if dist < rcut * (1 + eps_ver)
                                Num(alpha) = Num(alpha) + 1;
                                Neighborlist(alpha, Num(alpha), 1) = beta;

                                for i = 1:dim
                                    Neighborlist(alpha, Num(alpha), i + 1) = n;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        %disp('neighborlist not being updated');
    end
end
