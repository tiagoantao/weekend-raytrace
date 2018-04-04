// Chapter 1 - output an image
const get_red_index = (height, width, x, y) => {
    return (height - y)*(width*4) + x*4
}


const hello_world_1 = (canvas_id) => {
    const canvas = document.getElementById(canvas_id)
    const ctx = canvas.getContext("2d")

    const image_data = ctx.getImageData(0, 0, canvas.width, canvas.height)
    const data = image_data.data
    for (var y=canvas.height-1; y>=0; y--) {
        for (var x=0; x<canvas.width; x++) {
            const red_index = get_red_index(canvas.height, canvas.width, x, y)
            data[red_index] = Math.round(255 * x / canvas.width)
            data[red_index + 1] = y; Math.round(255 * y / canvas.height)
            data[red_index + 2] = 100
            data[red_index + 3] = 255
        }
    }
    ctx.putImageData(image_data, 0, 0)
}

// Chapter 2 - vectors (3 components)


const vec_read = (in_s) => {
}


const vec_write = (v, out) => {
}


const vec_squared_length = (v) => v.reduce((acu, cur) => acu + cur*cur, 0)


const vec_length = (v) => Math.sqrt(vec_squared_length(v))


const vec_unit = (v) => {
    const l = vec_length(v)
    return v.map(c => c/l)
}


const vec_sum = (v1, v2) => v1.map((v, i) => v + v2[i])


const vec_sub = (v1, v2) => v1.map((v, i) => v - v2[i])


const vec_mul = (v1, v2) => v1.map((v, i) => v*v2[i])


const vec_div = (v1, v2) => v1.map((v, i) => v/v2[i])


const vec_mul_num = (v, n) => v.map(val => val*n)


const vec_div_num = (v, n) => v.map(val => val/n)


const vec_dot = (v1, v2) => v1.reduce((acu, cur, i) => acu + cur*v2[i], 0)


const vec_cross = (v1, v2) => [
    v1[1]*v2[2] - v1[2]*v2[1],
    -(v1[0]*v2[2] - v1[2]*v2[0]),
    v1[0]*v2[1] - v1[1]*v2[0],
]


// Chapter 3 - Rays, a simple camera, and background


const point_at_parameter = (t, ray) => vec_sum(
    ray.origin,
    vec_mul_num(ray.direction, t)
)


const color_3 = (ray) => {
    const unit_direction = vec_unit(ray.direction)
    const t = 0.5*(unit_direction[1] + 1.0)
    const tm1 = 1.0 - t
    return vec_sum(
        [tm1, tm1, tm1],
        vec_mul_num([0.5, 0.7, 1.0], t)
    )
}


const hello_world_3 = (canvas_id, my_color=color_3) => {
    const canvas = document.getElementById(canvas_id)
    const ctx = canvas.getContext("2d")

    const width = canvas.width
    const height = canvas.height
    const image_data = ctx.getImageData(0, 0, width, height)
    const data = image_data.data

    const lower_left_corner = [-2.0, -1.0, -1.0]
    const horizontal = [4.0, 0.0, 0.0]
    const vertical = [0.0, 2.0, 0.0]
    const origin = [0.0, 0.0, 0.0]
    for (var y=height-1; y>=0; y--) {
        for (var x=0; x<width; x++) {
            const u = x / width
            const v = y / height
            const direction = vec_sum(
                lower_left_corner,
                vec_sum(
                    vec_mul_num(horizontal, u),
                    vec_mul_num(vertical, v)
                ))
            const ray = {origin, direction}
            const color = my_color(ray)
            const red_index = get_red_index(height, width, x, y)
            data[red_index] = Math.round(259.9*color[0])
            data[red_index + 1] = Math.round(259.9*color[1])
            data[red_index + 2] = Math.round(259.9*color[2])
            data[red_index + 3] = 255
        }
    }
    ctx.putImageData(image_data, 0, 0)
}


// Chapter 4 - Adding a sphere

const z_sphere = {
    center: [0, 0, -1],
    radius: 0.5
}

const compute_discriminant = (sphere, ray) => {
    const oc = vec_sub(ray.origin, sphere.center)
    const a = vec_dot(ray.direction, ray.direction)
    const b = 2.0 * vec_dot(oc, ray.direction)
    const c = vec_dot(oc, oc) - sphere.radius*sphere.radius
    const discriminant = b*b - 4*a*c
    return {a, b, c, discriminant}
}


const hit_sphere_4 = (sphere, ray) => {
    return compute_discriminant(sphere, ray).discriminant > 0
}


const color_4 = (ray) => {
    if (hit_sphere_4(z_sphere, ray)) {
        return [1, 0, 0]
    }
    return color_3(ray)
}


// Chapter 5 - Part 1 Surface Normals


const hit_sphere_5 = (sphere, ray) => {
    const discriminant = compute_discriminant(sphere, ray)
    return discriminant.discriminant < 0 ?
        -1.0
        :
        (-discriminant.b - Math.sqrt(discriminant.discriminant)) / (
            2.0*discriminant.a)
}


const color_5_1 = (ray) => {
    const t1 = hit_sphere_5(z_sphere, ray)
    if (t1 > 0.0) {
        const N = vec_unit(
            vec_sub(
                point_at_parameter(t1, ray),
                [0, 0, -1]))
        return vec_mul_num([N[0]+1, N[1]+1, N[2]+1], 0.5)
    }
    /*
    const unit_direction = vec_unit(ray.direction)
    const t2 = 0.5 * (unit_direction[1] + 1.0)
    return vec_sum(
        vec_mul_num([1, 1, 1], (1 - t2)),
        vec_mul_num([0.5, 0.7, 1.0], t2))
    */
    return color_3(ray)
}


// Chapter 5 - Part 2 Multiple objects

//Poor human's polymorphism
const hit_poly = (obj, ray, tmin, tmax) => {
    if (obj instanceof Array) {
        return hit_list(obj, ray, tmin, tmax)
    }
    return hit_sphere(obj, ray, tmin, tmax)
}

const hit_sphere = (sphere, ray, tmin, tmax) => {
    const d = compute_discriminant(sphere, ray)
    if (d <= 0) return false
    const temp1 = (-d.b - Math.sqrt(d.b*d.b - 4*d.a*d.c)) / (2*d.a)
    if (temp1 < tmax && temp1 > tmin) {
        const p = point_at_parameter(temp1, ray)
        return {
            t: temp1,
            p,
            normal: vec_div_num(
                vec_sub(p, sphere.center),
                sphere.radius
            ),
            mat: sphere.mat  // Chapter 8... 
        }
    }

    const temp2 = (-d.b + Math.sqrt(d.b*d.b - 4*d.a*d.c)) / (2*d.a)
    if (temp2 < tmax && temp2 > tmin) {
        const p = point_at_parameter(temp2, ray)
        return {
            t: temp2,
            p,
            normal: vec_div_num(
                vec_sub(p, sphere.center),
                sphere.radius
            ),
            mat: sphere.mat  // Chapter 8... 
        }
    }
    
    return false
}


const hit_list = (list, ray, tmin, tmax) => {
    let temp_rec = false
    let closest_so_far = tmax
    for (let elem of list) {
        const erec = hit_poly(elem, ray, tmin, closest_so_far)
        if (erec) {
            closest_so_far = erec.t
            temp_rec = erec
        }
    }
    return temp_rec
}

const color_5_2 = (ray, world) => {
    const rec = hit_list(world, ray, 0.0, Number.MAX_VALUE)
    return rec?
        vec_mul_num(
            [rec.normal[0]+1, rec.normal[1]+1, rec.normal[2]+1],
            0.5)
        :
        color_3(ray)
}


// Chapter 6 - Antialiasing


const simple_camera = {
    lower_left_corner: [-2.0, -1.0, -1.0],
    horizontal: [4.0, 0.0, 0.0],
    vertical: [0.0, 2.0, 0.0],
    origin:  [0.0, 0.0, 0.0]
}

const get_ray = (camera, u, v) => {
    const direction = vec_sum(
        camera.lower_left_corner,
        vec_sum(
            vec_mul_num(camera.horizontal, u),
            vec_mul_num(camera.vertical, v)
        ))
    const ray = {origin: camera.origin, direction}
    return ray
}

const hello_world_6 = (canvas_id, camera, num_samples, my_color) => {
    const canvas = document.getElementById(canvas_id)
    const ctx = canvas.getContext("2d")

    const width = canvas.width
    const height = canvas.height
    const image_data = ctx.getImageData(0, 0, width, height)
    const data = image_data.data

    for (var y=height-1; y>=0; y--) {
        for (var x=0; x<width; x++) {
            let color = [0, 0, 0]
            for (let sample=0; sample<num_samples; sample++) {
                const u = (x + Math.random()) / width
                const v = (y + Math.random()) / height
            
                const ray = get_ray(camera, u, v)
                color = vec_sum(my_color(ray), color)

            }
            color = vec_div_num(color, num_samples)
            const red_index = get_red_index(height, width, x, y)
            data[red_index] = Math.round(259.9*color[0])
            data[red_index + 1] = Math.round(259.9*color[1])
            data[red_index + 2] = Math.round(259.9*color[2])
            data[red_index + 3] = 255
        }
    }
    ctx.putImageData(image_data, 0, 0)
}


// Chapter 7 - Diffuse materials


const random_in_unit_sphere = () => {
    let v = [0, 0, 0]
    do {
        v = [2*Math.random() - 1, 2*Math.random() - 1, 2*Math.random() - 1]
    } while (vec_squared_length(v) >= 1)
    return v
}


const color_7 = (ray, world, acne=true) => {
    const rec = hit_list(world, ray, acne? 0.0: 0.01, Number.MAX_VALUE)
    if (rec) {
        const target = vec_sum(
            vec_sum(rec.p, rec.normal),
            random_in_unit_sphere()
        )
        return vec_mul_num(
            color_7({origin: rec.p, direction: vec_sub(target, rec.p)},
                    world),
            0.5)
    }
    else {
        return color_3(ray)
    }
}


const color_7_gamma = (ray, world) => color_7(ray, world).map(c => Math.sqrt(c))


// Chapter 8 - Metal


const scatter_lambertian = (albedo, ray_in, rec) => {
    const scat = [0, 0, 0]
    const target = vec_sum(
        vec_sum(rec.p, rec.normal),
        random_in_unit_sphere())
    const scatter = {origin: rec.p, direction: vec_sub(target, rec.p)}
    const attenuation = albedo
    return {attenuation, scatter}
}


const reflect = (v, n) => vec_sub(v, vec_mul_num(n, 2*vec_dot(v, n)))

const scatter_metal = (albedo, fuz, ray_in, rec) => {
    const reflected = reflect(vec_unit(ray_in.direction), rec.normal)
    const scatter = {
        origin: rec.p,
        direction: vec_sum(reflected,
                           vec_mul_num(random_in_unit_sphere(), fuz))}
    const attenuation = albedo
    return vec_dot(scatter.direction, rec.normal) > 0 ?
        {attenuation, scatter} : false
}


const scatter_poly = (mat, ray_in, rec) => {
    if (mat.mat === 'metal') {
        return scatter_metal(mat.albedo, mat.fuz, ray_in, rec)
    }
    else if (mat.mat === 'dielectric') {
        return scatter_dielectric(mat.ri, ray_in, rec)
    }
    return scatter_lambertian(mat.albedo, ray_in, rec)
}

const color_8 = (ray, world, depth=0) => {
    const rec = hit_list(world, ray, 0.001, Number.MAX_VALUE)
    if (rec) {
        const scatter = scatter_poly(rec.mat, ray, rec)
        if (scatter && depth < 50) { //50 on book
            return vec_mul(
                color_8(scatter.scatter, world, depth + 1),
                scatter.attenuation)
        }
        else {
            return [0, 0, 0]
        }
    }
    else {
        return color_3(ray)
    }
}


// Chapter 9 - Dielectrics

const refract = (v, n, ni_over_nt) => {
    const uv = vec_unit(v)
    const dt = vec_dot(uv, n)
    const discriminant = 1 - ni_over_nt*ni_over_nt*(1 - dt*dt)
    //console.log(ni_over_nt, uv, dt, n, discriminant)
    return discriminant > 0 ?
        vec_sub(
            vec_mul_num(
                vec_sub(uv, vec_mul_num(n, dt)), ni_over_nt),
            vec_mul_num(n, Math.sqrt(discriminant)))
        :
        false
}


const schlick = (cosine, ri) => {
    const r0 = (1-ri) / (1+ri)
    const r02 = r0*r0
    return r02 * (1-r02)*Math.pow((1 - cosine), 5)
}


const scatter_dielectric = (ri, ray_in, rec) => {
    let outward_normal
    let ni_over_nt
    const attenuation = [1, 1, 1]
    let cosine
    const reflected = reflect(ray_in.direction, rec.normal)

    if (vec_dot(ray_in.direction, rec.normal) > 0) {
        outward_normal = rec.normal.map(c => -c)
        ni_over_nt = ri
        cosine = ri * vec_dot(ray_in.direction, rec.normal) / vec_length(ray_in.direction)
    }
    else {
        outward_normal = rec.normal
        ni_over_nt = 1.0 / ri
        cosine = -vec_dot(ray_in.direction, rec.normal) / vec_length(ray_in.direction)
    }
    const refracted = refract(ray_in.direction, outward_normal, ni_over_nt)
    //console.log(refracted)
    const reflect_prob = refracted ?
          schlick(cosine, ri)
          :
          1
    //console.log(refracted, reflect_prob)
    return Math.random() < reflect_prob ?
        {attenuation, scatter: {origin: rec.p, direction: reflected}}
        :
        {attenuation, scatter: {origin: rec.p, direction: refracted}}
}


// Chapter 10 - Positionable camera

const standard_camera = (vfov, aspect) => {
    const theta = vfov * Math.PI / 180
    const half_height = Math.tan(theta/2)
    const half_width = aspect * half_height
    return {
        lower_left_corner: [-half_width, -half_height, -1.0],
        horizontal: [2*half_width, 0.0, 0.0],
        vertical: [0.0, 2*half_height, 0.0],
        origin:  [0.0, 0.0, 0.0]
    }
}

// boot


const world = [
    {center: [0, 0, -1], radius: 0.5},
    {center: [0, -100.5, -1], radius: 100}
]

const world_mat = [
    {center: [0, 0, -1], radius: 0.5,
     mat: {mat: 'lambertian', albedo: [0.8, 0.3, 0.3]}},
    {center: [0, -100.5, -1], radius: 100,
     mat: {mat: 'lambertian', albedo: [0.8, 0.8, 0]}},
    {center: [1, 0, -1], radius: 0.5,
     mat: {mat: 'metal', albedo: [0.8, 0.6, 0.2], fuz: 0}},
    {center: [-1, 0, -1], radius: 0.5,
     mat: {mat: 'metal', albedo: [0.8, 0.8, 0.8], fuz: 0}}
]


const world_dielectric = [
    {center: [0, 0, -1], radius: 0.5,
     mat: {mat: 'lambertian', albedo: [0.1, 0.2, 0.5]}},
    {center: [0, -100.5, -1], radius: 100,
     mat: {mat: 'lambertian', albedo: [0.8, 0.8, 0]}},
    {center: [1, 0, -1], radius: 0.5,
     mat: {mat: 'metal', albedo: [0.8, 0.6, 0.2], fuz: 0}},
    {center: [-1, 0, -1], radius: 0.5,
    mat: {mat: 'dielectric', ri: 1.5}},
    {center: [-1, 0, -1], radius: -0.45,
     mat: {mat: 'dielectric', ri: 1.5}}
]

const _R = Math.cos(Math.PI / 4)
const world_camera_1 = [
    {center: [-_R, 0, -1], radius: _R,
     mat: {mat: 'lambertian', albedo: [0, 0, 01]}},
    {center: [_R, 0, -1], radius: _R,
     mat: {mat: 'lambertian', albedo: [1, 0, 0]}}
]


const start = () => {
    //hello_world_3('ray1', color_3)
    //hello_world_3('ray1', color_4)
    //hello_world_3('ray1', color_5_1)
    //hello_world_3('ray1', (ray) => color_5_2(ray, world))


    /* antialias in book is 100, not 10 */
    
    //hello_world_6('ray1', simple_camera, 1, (ray) => color_5_2(ray, world))
    //hello_world_6('ray2', simple_camera, 10, (ray) => color_5_2(ray, world))

    //hello_world_6('ray1', simple_camera, 1, (ray) => color_7(ray, world))
    //hello_world_6('ray2', simple_camera, 1, (ray) => color_7_gamma(ray, world))
    //hello_world_6('ray3', simple_camera, 1, (ray) => color_7_gamma(ray, world, false))
    //hello_world_6('ray4', simple_camera, 10, (ray) => color_7_gamma(ray, world, false))

    //hello_world_6('ray1', simple_camera, 5, (ray) => color_8(ray, world_mat))

    //hello_world_6('ray1', simple_camera, 20, (ray) => color_8(ray, world_dielectric))

    const canvas_camera_1 = 'ray1'
    const canvas = document.getElementById(canvas_camera_1)
    const ctx = canvas.getContext("2d")
    const width = canvas.width
    const height = canvas.height
    hello_world_6('ray1', standard_camera(90, width/height), 20, (ray) => color_8(ray, world_camera_1))
    
}
