//todo - need to pull out localStorage service so not internal so can be tested

describe('clotho.core Collector tests', function() {
	var factory;

	// excuted before each "it()" is run.
	beforeEach(function() {
		// load the module
		module('clotho.core');

		// inject your factory for testing
		inject(function(Collector) {
			factory = Collector;
		});

		//this mocks localStorage
		var store = {
			123456789: {},
			asdfghjkl: {},
			qwertyuiop: {}
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

	it('should have a storeModel function', function() {
		expect(angular.isFunction(factory.storeModel)).toBe(true);
	});

	it('should have a retrieveModel function', function() {
		expect(angular.isFunction(factory.retrieveModel)).toBe(true);
	});

	it('should should return the whole collector object without arguments', function() {
		expect(angular.isObject(factory.collector)).toBe(true);
	});

	it('should return items by key when call retrieveModel', function() {
		var result = factory.retrieveModel(123456789);

		expect(result.length == 1 && angular.isObject(result)).toBe(true);
	});

	it('should add an item when calling storeModel and retrieve it with retrieveModel', function() {
		var myObj = {hi: 'there'};
		factory.storeModel('12849124219', myObj);

		expect(factory.retrieveModel('12849124219')).toEqual(myObj);
	});
});