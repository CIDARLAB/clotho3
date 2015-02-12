describe("query", function(){
    "use strict";
    //create test objects so tests are stable
    beforeEach(function(done){
        var thiz = this;
        Clotho.login("read", "read")
        .then(function(success){
            var testitems = [];
            //superschema tester
            testitems.push({id:"org.clothocad.Tester",
                            schema:"org.clothocad.core.schema.ClothoSchema"});
            //subschema specifictester
            testitems.push({id:"org.clothocad.SpecificTester",
                            schema:"org.clothocad.core.schema.ClothoSchema",
                            superClass:"org.clothocad.Tester"});
            //one tester, one specifictester
            testitems.push({name:"A Tester", 
                            schema:"org.clothocad.Tester"});
            testitems.push({name:"A Tester", 
                            schema:"org.clothocad.SpecificTester"});

            return Clotho.create(testitems);
        })
        .then(function(data){
            thiz.ids = data;
        })
        .finally(function(){done();});
    });

    afterEach(function(done) {
        Clotho.destroy(this.ids)
        .then(function(success){return Clotho.logout();})
        .finally(function(){done();});
    });

    it("finds objects that satisfy the description", function(done) {
        Clotho.query({name:"A Tester"})
        .then(function(data) {
            expect(data.length).toEqual(2);
            for (var i = 0; i < data.length; i++){
                expect(data[i].name).toEqual("A Tester");
            }
        })
        .finally(function(){done();});
    });
        
    it("finds subtypes of schemas", function(done) {
        Clotho.query({schema:"org.clothocad.Tester"})
        .then(function(data) {
            expect(data.length).toEqual(2);
        })
        .finally(function(){done();});
    });

    it("does not find supertypes of schemas", function(done) {
        Clotho.query({schema:"org.clothocad.SpecificTester"})
        .then(function(data) {
            expect(data.length).toEqual(1);
        })
        .finally(function(){done();});
    });

    it("returns empty list if 0 results", function(done) {
        Clotho.query({schema:"org.clothocad.SirNotAppearingInThisTest"})
        .then(function(data) {
            expect(data.length).toEqual(0);
        })
        .finally(function(){done();});
    });

    describe("queryOne", function(){
        it("returns one result", function() {
           
            Clotho.queryOne({schema:"org.clothocad.Tester"})
            .then(function(data) {
                expect(data.schema).toEqual("org.clothocad.Tester");
            })
            .finally(function(){done();});
            
        });

        it("returns nothing if no results", function() {
            Clotho.queryOne({schema:"org.clothocad.SirNotAppearingInThisTest"})
            .then(function(data) {
                expect(data).toBeUndefined();
            })
            .finally(function(){done();});
        }); 
    }); 
});
