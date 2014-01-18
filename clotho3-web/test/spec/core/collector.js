//todo - need to pull out localStorage service so not internal so can be tested

describe('clotho.core Collector', function() {
	var factory,
		window;

	// excuted before each "it()" is run.
	beforeEach(function() {
		// load the module
		module('clotho.core');

		// inject your factory for testing
		inject(function(Collector, $window) {
			factory = Collector;
			window = $window
		});

		//this mocks localStorage
		var store = {
			"clotho_123456789": {"hey" : "there"},
			"clotho_qwertyuiop": {"nothing" : "here"}
		};

		spyOn(localStorage, 'getItem').andCallFake(function(key) {
			return store[key];
		});

		spyOn(localStorage, 'setItem').andCallFake(function(key, value) {
			return store[key] = value + '';
		});

		spyOn(localStorage, 'clear').andCallFake(function() {
			store = {};
		});

		spyOn(Object, 'keys').andCallFake(function(value) {
			var keys=[];

			for(var key in store) {
				keys.push(key);
			}

			return keys;
		});
	});


	it('should should return the whole collector object without arguments', function() {
		expect(angular.isObject(factory.collector)).toBe(true);
	});


	describe('#retrieveModel', function() {

		it('should have a retrieveModel function', function() {
			expect(angular.isFunction(factory.retrieveModel)).toBe(true);
		});

		it('should return items by key when call retrieveModel', function() {
			var result = factory.retrieveModel('123456789');

			expect(angular.isObject(result)).toBe(true);
		});

		it('should prefix a key with clotho_', function() {
			var direct = window.localStorage.getItem('clotho_qwertyuiop');
			var indirect = factory.retrieveModel('qwertyuiop');
			expect(indirect).toEqual(direct);
		});

	});

	describe('#storeModel', function() {

		it('should have a storeModel function', function() {
			expect(angular.isFunction(factory.storeModel)).toBe(true);
		});

		it('should add an item when calling storeModel and retrieve it with retrieveModel', function() {
			var myObj = {hi: 'there'};
			factory.storeModel('12849124219', myObj);

			expect(factory.retrieveModel('12849124219')).toEqual(myObj);
		});

	});
});