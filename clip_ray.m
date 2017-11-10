function [r_0c, r_1c] = clip_ray(r_0, r_1, a)
% Function clip_ray: Clips a ray. At this point, I've rewritten the
% function 4 times to gracefully handle pathological cases. Please bear
% with me if the documentation seems lazy, I'd highly recommend that you
% consider the situation I'm in.



% Before we begin, let's express the parametric equation in terms of t and
% v.
v = 1 / norm(r_1 - r_0) * (r_1 - r_0);

if max(abs(r_0)) > a
    
    min_t = Inf;
    
    % Check intersections for each axis...
    for i = 1:3
        
        % Only worth doing if the ray component will intersect at all. So
        % we check if the v component will be zero...
        
        if v(i) ~= 0
            pos_int_t = ( a - r_0(i)) / v(i);
            neg_int_t = (-a - r_0(i)) / v(i);
            
            pos_point = r_0 + pos_int_t * v;
            neg_point = r_0 + neg_int_t * v;
            
            % These will be used to check which POI is closest on the box,
            % if they hit it at all...
            t_cmp_1 = Inf;
            t_cmp_2 = Inf;
            
            % Not useful if they're not in the box domain...
            if round(max(abs(pos_point)) - 1e-10) <= a && pos_int_t >= 0
                t_cmp_1 = pos_int_t;
            end
            
            if round(max(abs(neg_point)) - 1e-10) <= a && neg_int_t >= 0
                t_cmp_2 = neg_int_t;
            end
            % The small number subtraction above will account for floating
            % point errors. It's still very tiny, so it's a safe
            % approximation. If it's anything, we added the condition that
            % the intersection times must be positive...
            
            % Now we check which side of the box is closest...
            min_t = min(min_t, min(t_cmp_1, t_cmp_2));
            
        end
            
            
            
    end
    
    r_0 = r_0 + min_t * v;
    
end

% Now we do it for the other side. Here, we simply flip v and take r_1 to
% be the starting point. Then we do the same thing. 
v_f = -v;

if max(abs(r_1)) > a
    
    min_t = Inf;
    
    % Check intersections for each axis...
    for i = 1:3
        
        % Only worth doing if the ray component will intersect at all. So
        % we check if the v component will be zero...
        
        if v_f(i) ~= 0
            pos_int_t = ( a - r_1(i)) / v_f(i);
            neg_int_t = (-a - r_1(i)) / v_f(i);
            
            pos_point = r_1 + pos_int_t * v_f;
            neg_point = r_1 + neg_int_t * v_f;
            
            % These will be used to check which POI is closest on the box,
            % if they hit it at all...
            t_cmp_1 = Inf;
            t_cmp_2 = Inf;
            
            % Not useful if they're not in the box domain...
            if round(max(abs(pos_point)) - 1e-10) <= a && pos_int_t >= 0
                t_cmp_1 = pos_int_t;
            end
            
            if round(max(abs(neg_point)) - 1e-10) <= a && neg_int_t >= 0
                t_cmp_2 = neg_int_t;
            end
            
            % Now we check which side of the box is closest...
            min_t = min(min_t, min(t_cmp_1, t_cmp_2));
            
        end
            
            
            
    end
    
    r_1 = r_1 + min_t * v_f;
    
end


r_0c = r_0;
r_1c = r_1;

end
