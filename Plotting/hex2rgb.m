function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1) == '#'
        hex = hex(2:end);
    end
    if length(hex) ~= 6
        error('Input must be a 7-character string including "#" or a 6-character string.');
    end
    r = hex2dec(hex(1:2));
    g = hex2dec(hex(3:4));
    b = hex2dec(hex(5:6));
    rgb = [r, g, b] / 255;
end
