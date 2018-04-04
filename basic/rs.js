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

//Poor human polymorphism
const poly_hit = (obj, ray, tmin, tmax) => {
    if (obj instanceof Array) {
        return list_hit(obj, ray, tmin, tmax)
    }
    return sphere_hit(obj, ray, tmin, tmax)
}

const sphere_hit = (sphere, ray, tmin, tmax) => {
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
            )
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
            )
        }
    }
    
    return false
}


const list_hit = (list, ray, tmin, tmax) => {
    let temp_rec = false
    let closest_so_far = tmax
    for (let elem of list) {
        const erec = poly_hit(elem, ray, tmin, closest_so_far)
        if (erec) {
            closest_so_far = erec.t
            temp_rec = erec
        }
    }
    return temp_rec
}

const color_5_2 = (ray, world) => {
    const rec = list_hit(world, ray, 0.0, Number.MAX_VALUE)
    return rec?
        vec_mul_num(
            [rec.normal[0]+1, rec.normal[1]+1, rec.normal[2]+1],
            0.5)
        :
        color_3(ray)
}


// Chapter 6 - Antialiasing


const default_camera = {
    lower_left_corner: [-2.0, -1.0, -1.0],
    horizontal: [4.0, 0.0, 0.0],
    vertical: [0.0, 2.0, 0.0],
    origin:  [0.0, 0.0, 0.0]
}

const camera_get_ray = (camera, u, v) => {
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
            
                const ray = camera_get_ray(camera, u, v)
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



// boot


const world = [
    {center: [0, 0, -1], radius: 0.5},
    {center: [0, -100.5, -1], radius: 100}
]

const start = () => {
    //hello_world_3('ray1', color_3)
    //hello_world_3('ray1', color_4)
    //hello_world_3('ray1', color_5_1)
    //hello_world_3('ray1', (ray) => color_5_2(ray, world))
    hello_world_6('ray1', default_camera, 1, (ray) => color_5_2(ray, world))
    hello_world_6('ray2', default_camera, 10, (ray) => color_5_2(ray, world))
}
