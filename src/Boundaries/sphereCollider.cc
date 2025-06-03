#include <vector>
#include "../Physics/physics.hh"
#include "../DataBase/nodeList.hh"
#include "collider.hh"

// Base class for Grid Boundaries
template <int dim>
class SphereCollider : public Collider<dim> {
protected:
    Lin::Vector<dim> position;
    double radius,elasticity;

    Lin::Vector<dim>
    normal(Lin::Vector<dim>& otherPosition) {
        return (otherPosition-position).normal();
    }
public:
    using Vector=Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    SphereCollider(Vector& position, double radius, double elasticity) : 
        position(position), radius(radius), elasticity(elasticity) {
        if (elasticity > 1.0)
            std::cerr << "elasticity > 1. was this intentional?" << std::endl;
    }

    SphereCollider(Vector& position, double radius) : 
        SphereCollider(position,radius,1.0) {}
    

    virtual ~SphereCollider() {}

    virtual void
    ApplyBoundaries(State<dim>* state, NodeList* nodeList) override {  
        int numNodes            = state->size();
        VectorField* positions  = nodeList->getField<Vector>("position");
        VectorField* velocities = nodeList->getField<Vector>("velocity");
        ScalarField* radii      = nodeList->getField<double>("radius");

        if (radii!= nullptr) {
            for (int i=0; i<numNodes; ++i) {
                Vector pos = positions->getValue(i);
                double rad = radii->getValue(i);
                if (Inside(pos,rad)) {
                    Vector v1 = velocities->getValue(i);
                    Vector n  = (pos-position).normal();
                    Vector v2 = v1 - 2.0*(v1*n)*n;
                    velocities->setValue(i,v2*elasticity);                 // reflect across the normal
                    Vector newPos = pos + (0.1*rad)*n;  // move it out of the boundary
                    positions->setValue(i,newPos);
                }
            }
            //physics->PushState(state);
            /* i'm changing boundaries to modify the nodeList directly for now rather than the state
               this may not be a good idea for full generality, but at the moment, there is just a lot
               of copying back and forth to and from the state and that incurs a cost. */
        }
        else
        {
            std::cerr << "radius not found in nodeList. was this intentional?" << std::endl;
        }
    }

    bool 
    Inside(Lin::Vector<dim>& otherPosition, double otherRadius) override {
        double distance = (position-otherPosition).magnitude();
        if (distance <= (radius + otherRadius))
            return true;
        else
            return false;
    }

};