
// Chapter 1 - output an image
const get_red_index = (height, width, x, y) => {
    return (height - y)*(width*4) + x*4
}


const hello_world_1 = () => {
    const canvas = document.getElementById("ray")
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


const point_at_parameter = (t, origin, direction) => vec_sum(
    origin,
    vec_mul_num(direction, t)
)


const color_3 = (origin, direction) => {
    const unit_direction = vec_unit(direction)
    const t = 0.5*(unit_direction[1] + 1.0)
    const tm1 = 1.0 - t
    return vec_sum(
        [tm1, tm1, tm1],
        vec_mul_num([0.5, 0.7, 1.0], t)
    )
}


const hello_world_3 = (my_color=color_3) => {
    const canvas = document.getElementById("ray")
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
            const color = my_color(origin, direction)
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

const compute_discriminant = (center, radius, ray_origin, ray_direction) => {
    const oc = vec_sub(ray_origin, center)
    const a = vec_dot(ray_direction, ray_direction)
    const b = 2.0 * vec_dot(oc, ray_direction)
    const c = vec_dot(oc, oc) - radius*radius
    const discriminant = b*b - 4*a*c
    return {a, b, c, discriminant}
}


const hit_sphere_4 = (center, radius, ray_origin, ray_direction) => {
    return compute_discriminant(
        center, radius, ray_origin, ray_direction).discriminant > 0
}


const color_4 = (origin, direction) => {
    if (hit_sphere_4([0, 0, -1], 0.5, origin, direction)) {
        return [1, 0, 0]
    }
    return color_3(origin, direction)
}


// Chapter 5 - Surface Normals and multiple objects


const hit_sphere_5 = (center, radius, ray_origin, ray_direction) => {
    const discriminant = compute_discriminant(
        center, radius, ray_origin, ray_direction)
    return discriminant.discriminant < 0 ?
        -1.0 :
        (-discriminant.b - Math.sqrt(discriminant.discriminant)) / (
            2.0*discriminant.a)
}


const color_5 = (origin, direction) => {
    const t1 = hit_sphere_5([0, 0, -1], 0.5, origin, direction)
    if (t1 > 0.0) {
        const N = vec_unit(
            vec_sub(
                point_at_parameter(t1, origin, direction),
                [0, 0, -1]))
        return vec_mul_num([N[0]+1, N[1]+1, N[2]+1], 0.5)
    }
    const unit_direction = vec_unit(direction)
    const t2 = 0.5 * (unit_direction[1] + 1.0)
    return vec_sum(
        vec_mul_num([1, 1, 1], (1 - t2)),
        vec_mul_num([0.5, 0.7, 1.0], t2))
}

// boot


const start = () => {
    hello_world_3(color_5)
}
